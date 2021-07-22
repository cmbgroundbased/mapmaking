//! My version of the distribution of processes within the threads. Nowdays, It is works only for single node

use std::{sync::{Arc, Mutex}, thread};
use std::sync::mpsc;

enum Message {
    NewJob(Job),
    Terminate,
}

pub struct ThreadPool {
    workers: Vec<Worker>,
    sender: mpsc::Sender<Message>,
}


type Job = Box<dyn FnOnce() + Send + 'static>;


impl ThreadPool {
    pub fn new(size: usize) -> Self {
        assert!(size > 0);

        let (tx, rx) = mpsc::channel();

        let rx = Arc::new(Mutex::new(rx));

        let mut workers = Vec::with_capacity(size);

        for id in 0..size{
            workers.push(Worker::new(id, Arc::clone(&rx)));
        }

        ThreadPool {
            workers,
            sender: tx,
        }
    }

    // we need Send to transfer the closure from one thread
    // to another and 'static because we don’t know how long
    // the thread will take to execute.
    // pub fn spawn<F, T>(f: F) -> JoinHandle<T>
    //     where F: FnOnce()  + Send + 'static, T: Send + 'static {

    // }
    pub fn execute<F>(&self, f: F) where F: FnOnce() + Send + 'static {
        let job = Box::new(f);
        self.sender.send(Message::NewJob(job)).unwrap();

    }
}

struct Worker {
    id: usize,
    thread: Option<thread::JoinHandle<()>>,
}

impl Worker {
    fn new(id: usize, rx: Arc<Mutex<mpsc::Receiver<Message>>>) -> Self {
        let thread = thread::spawn(move || 
                loop {
                    let message = rx.lock().unwrap().recv().unwrap();
                    match message {
                        Message::NewJob(job) => {

                            //println!("Worker {} got a job; executing.", id);
                            job();

                        }

                        Message::Terminate => {

                            //println!("Worker {} was told to terminate.", id);
                            break;

                        }
                    }
                });

        Worker {
            id,
            thread: Some(thread),
        }
    }
}

impl Drop for ThreadPool {
    fn drop(&mut self) {
        //println!("Sending terminate message to all workers.");

        for _ in &self.workers {
            self.sender.send(Message::Terminate).unwrap();
        }

        //println!("Shutting down all workers.");

        for worker in &mut self.workers {
            //println!("Shutting down worker {}", worker.id);
            if let Some(t) = worker.thread.take() {
                t.join().unwrap();
            }
        }
    }
}
