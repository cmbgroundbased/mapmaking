extern crate num_cpus;

use std::fmt;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{channel, Receiver, Sender};
use std::sync::{Arc, Condvar, Mutex};
use std::thread;

// enum Message {
//     NewJob(Job<'static>),
//     Terminate,
// }

pub struct ThreadPool {
    sender: Sender<Job<'static>>,
    workers: Arc<Worker>,
}

struct Sentinel<'a> {
    shared_data: &'a Arc<Worker>,
    active: bool,
}

#[derive(Clone, Default)]
pub struct Builder {
    num_threads: Option<usize>,
    thread_name: Option<String>,
    thread_stack_size: Option<usize>,
}

struct Worker {
    name: Option<String>,
    job_receiver: Mutex<Receiver<Job<'static>>>,
    empty_trigger: Mutex<()>,
    empty_condvar: Condvar,
    join_generation: AtomicUsize,
    queued_count: AtomicUsize,
    active_count: AtomicUsize,
    max_thread_count: AtomicUsize,
    panic_count: AtomicUsize,
    stack_size: Option<usize>,
}


trait FnBox {
    fn call_box(self: Box<Self>);
}

impl<F: FnOnce()> FnBox for F {
    fn call_box(self: Box<F>) {
        (*self)();
    }
}

type Job<'a> = Box<dyn FnBox + Send + 'a>;

/******* The implementation of the Sentinel's metods */

impl<'a> Sentinel<'a> {
    fn new(shared_data: &'a Arc<Worker>) -> Sentinel<'a> {
        Sentinel {
            shared_data,
            active: true,
        }
    }

    fn cancel(mut self){
        self.active = false;
    }
}

impl<'a> Drop for Sentinel<'a> {
    fn drop(&mut self) {
        if self.active {
            self.shared_data.active_count.fetch_sub(1, Ordering::SeqCst);
            if thread::panicking() {
                self.shared_data.panic_count.fetch_add(1, Ordering::SeqCst);
            }
            self.shared_data.no_work_nofity_all();
            spawn_in_pool(self.shared_data.clone())
        }
    }
}


/******** The implementation of the Builder */

impl Builder {
    pub fn new() -> Self {
        Builder {
            num_threads: None,
            thread_name: None,
            thread_stack_size: None,
        }
    }

    pub fn num_threads(mut self, num_threads: usize) -> Self {
        assert!(num_threads > 0);
        self.num_threads = Some(num_threads);
        self
    }

    pub fn thread_name(mut self, name: String) -> Self {
        self.thread_name = Some(name);
        self
    }

    pub fn thread_stack_size(mut self, size: usize) -> Self {
        self.thread_stack_size = Some(size);
        self
    }

    pub fn build(self) -> ThreadPool {
        let (tx, rx) = channel::<Job<'static>>();
        let num_threads = self.num_threads.unwrap_or_else(num_cpus::get);
        let shared_data = Arc::new( Worker {
            name: self.thread_name,
            job_receiver: Mutex::new(rx),
            empty_condvar: Condvar::new(),
            empty_trigger: Mutex::new(()),
            join_generation: AtomicUsize::new(0),
            queued_count: AtomicUsize::new(0),
            active_count: AtomicUsize::new(0),
            max_thread_count: AtomicUsize::new(num_threads),
            panic_count: AtomicUsize::new(0),
            stack_size: self.thread_stack_size,
        });

        for _ in 0..num_threads {
            spawn_in_pool(shared_data.clone());
        }

        ThreadPool {
            sender: tx,
            workers: shared_data,
        }
    }
}

/************************ The implementation of the Worker */

impl Worker {
    fn has_work(&self) -> bool {
        self.queued_count.load(Ordering::SeqCst) > 0 || self.active_count.load(Ordering::SeqCst) > 0
    }

    fn no_work_nofity_all(&self) {
        if !self.has_work() {
            *self
                .empty_trigger
                .lock()
                .expect("Unable to notify all joining threads");
            self.empty_condvar.notify_all();
        }
    }
}

/**************************The implementation of the ThreadPool */

impl ThreadPool {
    pub fn new(num_threads: usize) -> Self {
        Builder::new().num_threads(num_threads).build()
    }

    pub fn with_name(name: String, num_threads: usize) -> Self {
        Builder::new()
                .num_threads(num_threads)
                .thread_name(name)
                .build()
    }

    pub fn execute<F>(&self, job: F) where F: FnOnce() + Send + 'static, {
        self.workers.queued_count.fetch_add(1, Ordering::SeqCst);
        self.sender.send(Box::new(job))
                   .expect("ThreadPool::execute unable to send job into the queue");
    }

    pub fn queued_count(&self) -> usize {
        self.workers.queued_count.load(Ordering::Relaxed)
    }

    pub fn active_count(&self) -> usize {
        self.workers.active_count.load(Ordering::SeqCst)
    }

    pub fn max_count(&self) -> usize {
        self.workers.max_thread_count.load(Ordering::Relaxed)
    }

    pub fn panic_count(&self) -> usize {
        self.workers.panic_count.load(Ordering::Relaxed)
    }

    pub fn set_num_threads(&mut self, num_threads: usize) {
        assert!(num_threads >= 1);
        let prev_num_threads = self
                    .workers
                    .max_thread_count
                    .swap(num_threads, Ordering::Relaxed);
        if let Some(num_spawn) = num_threads.checked_sub(prev_num_threads) {
            for _ in 0..num_spawn {
                spawn_in_pool(self.workers.clone());
            }
        }
    }

    pub fn join(&self) {
        if self.workers.has_work() == false {
            return ();
        }

        let generation = self.workers.join_generation.load(Ordering::SeqCst);
        let mut lock = self.workers.empty_trigger.lock().unwrap();

        while generation == self.workers.join_generation.load(Ordering::Relaxed) && self.workers.has_work() {
            lock = self.workers.empty_condvar.wait(lock).unwrap();
        }

        // self.workers.join_generation.compare_and_swap(
        //     generation, generation.wrapping_add(1),
        //     Ordering::SeqCst,
        // );
        self.workers.join_generation.compare_exchange(
            generation,
            generation.wrapping_add(1),
            Ordering::SeqCst, Ordering::Relaxed).unwrap();
    }
}

impl Clone for ThreadPool {
    fn clone(&self) -> Self {
        ThreadPool {
            workers: self.workers.clone(),
            sender: self.sender.clone(),
        }
    }
}

impl Default for ThreadPool {
    fn default() -> Self {
        ThreadPool::new(num_cpus::get())
    }
}

impl fmt::Debug for ThreadPool {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("ThreadPool")
            .field("name", &self.workers.name)
            .field("queued_count", &self.queued_count())
            .field("active_count", &self.active_count())
            .field("max_count", &self.max_count())
            .finish()
    }
}

impl PartialEq for ThreadPool {
    fn eq(&self, other: &ThreadPool) -> bool {
        //let a: &Worker = &*self.workers;
        //let b: &Worker = &*other.workers;


        Arc::ptr_eq(&self.workers, &other.workers)
    }
}

impl Eq for ThreadPool {}

fn spawn_in_pool(worker: Arc<Worker>) {
    let mut builder = thread::Builder::new();

    if let Some(ref name) = worker.name {
        builder = builder.name(name.clone());
    }

    if let Some(ref stack_size) = worker.stack_size {
        builder = builder.stack_size(stack_size.to_owned());
    }

    builder.spawn(move || {
        let sentinel = Sentinel::new(&worker);

        loop {

            let thread_counter_val = worker.active_count.load(Ordering::Acquire);
            let max_thread_count_val = worker.max_thread_count.load(Ordering::Relaxed);
            if thread_counter_val >= max_thread_count_val {
                break;
            }
            let message = {
                let lock = worker
                    .job_receiver
                    .lock()
                    .expect("Worker thread unable to lock job_receiver");
                lock.recv()
            };

            let job = match message {
                Ok(job) => job,
                Err(..) => break,
            };

            worker.active_count.fetch_add(1, Ordering::SeqCst);
            worker.queued_count.fetch_sub(1, Ordering::SeqCst);

            job.call_box();

            worker.active_count.fetch_add(1, Ordering::SeqCst);
            worker.no_work_nofity_all();

        }

        sentinel.cancel();
    }).unwrap();
}
