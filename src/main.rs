extern crate mapmaking;
// pub mod threadpool;
pub mod threadpool2;
//use colored::Colorize;
use threadpool2::ThreadPool;

//use std::{slice::SliceIndex, sync::{Arc, Barrier, atomic::AtomicUsize, mpsc::channel}, thread, time::Duration};
use clap::{App, Arg};
mod directory;
use directory::DirStruct;
use std::path::Path;
use mapmaking::{Obs, sky};

use std::thread::sleep;
use gnuplot::{Figure, Caption, Color};
// use mapmaking::iteratorscustom::FloatIterator;

fn main() {
	let program = App::new("Strip MapMaker")
					 .version("0.1.0")
					 .author("Stefano Mandelli")
					 .about("MapMaking program to mitigate the systematic effects and recontruct the maps in Intensity and Polarization")
					 .arg(Arg::with_name("tod_path").short("t").long("tod_path").takes_value(true).help("Path of the TOD's direcotry"))
					 .arg(Arg::with_name("mc_id").short("i").long("mc_id").takes_value(true).help("Returns the map at the `mc_id` Monte Carlo iteration"))
					 .get_matches();

	const NUM_THREADS: usize = 1;
	const NUM_MC_ITER: usize = 50;
	let my_pool = ThreadPool::new(NUM_THREADS);


	let mc_id_set = ["43000000.0", "43000001.0", "43000002.0", "43000003.0", "43000004.0", "43000005.0", "43000006.0", "43000007.0", "43000008.0", "43000009.0", "43000010.0", "43000011.0",
                              "43000012.0", "43000013.0", "43000014.0", "43000015.0", "43000016.0", "43000017.0", "43000018.0", "43000019.0", "43000020.0", "43000021.0", "43000022.0", "43000023.0",
                              "43000024.0", "43000025.0", "43000026.0", "43000027.0", "43000028.0", "43000029.0", "43000030.0", "43000031.0", "43000032.0", "43000033.0", "43000034.0", "43000035.0",
                              "43000036.0", "43000037.0", "43000038.0", "43000039.0", "43000040.0", "43000041.0", "43000042.0", "43000043.0", "43000044.0", "43000045.0", "43000046.0", "43000047.0",
                              "43000048.0", "43000049.0"];

	let directory_path = program.value_of("tod_path").unwrap();
	//let (tx, rx) = channel(); (possiamo far tornare cose al thread 0, in modo tale da tener traccia di quello che sta iniziando/finendo)

	for _th in 0..NUM_MC_ITER {
		let path = String::from(directory_path.clone());
		let id = mc_id_set[_th];

		my_pool.execute(move || {
			split_mc(path, id);
			sleep(std::time::Duration::from_millis(500));
		});
	}
}

pub fn split_mc(p: String, id: &str) -> &str {

	let directory_tree = DirStruct::new(Path::new(&*p), String::from("1")).unwrap();
	let my_sky = sky::Sky::new();
	let t_sky = my_sky.get_t_sky();


	let _my_obs = directory_tree.create_observations(id, t_sky);
  
	// _my_obs.binning();
	// _my_obs.dummy_denoise();
	_my_obs.gls_denoise(1E-5, 10, 128);
	
	// _my_obs.atm_mitigation(1, 10, 1E-10, 128);
	id

}
