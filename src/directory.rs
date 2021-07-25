//! Create directory structure
use crate::Obs;
use std::fs::read_dir;
use std::path::{Path, PathBuf};
use colored::*;
use std::io;
use ndarray::Array1;
use ndarray_npy::NpzReader;
use std::fs::File;
use std::marker::PhantomData;

pub struct DirStruct <'a> {
    _outer_directory: Vec<PathBuf>,
    mc_path: Vec<PathBuf>,
    _mc_samples: usize,
    _start_date: String,
    _observation_time: String,
    entry_point: String,

    _phantom: &'a PhantomData<f32>,
}

impl <'a> DirStruct <'a> {
    pub fn new(path: &Path, obs_time: String) -> Result<Self, Box<dyn std::error::Error>> {
        let mut _mc_samples = 0;
        let mut _start_date = "";
        let mut _stop_date = "";
        let mut mc_vec: Vec<PathBuf> = Vec::new();
        let mut _outer_directory: Vec<PathBuf> = Vec::new();

        let entries = read_dir(path)?
            .map(|res| res.map(|e| e.path()))
            .collect::<Result<Vec<_>, io::Error>>()?;

        let mut sort = entries.clone();
        sort.sort();

        _start_date = match sort[0].to_str() {
            Some(v) => v,
            None => "Empty value",
        };

        let last_elem = sort.len();

        _stop_date = match sort[last_elem - 1].to_str() {
            Some(v) => v,
            None => "Empty value",
        };

        for mc_i in entries.clone() {
            // mc_i = prova<data>, prova<data2>, etc...
            _outer_directory.push(mc_i.clone());
            let mc = read_dir(mc_i)?
                .map(|res| res.map(|e| e.path()))
                .collect::<Result<Vec<_>, io::Error>>()?;
            // mc ci sono tutte le directory montecarlo
            _mc_samples = mc.len();
            for i in mc {
                if i.is_dir() {
                    let freq_mc = read_dir(i)?
                        .map(|res| res.map(|a| a.path()))
                        .collect::<Result<Vec<_>, io::Error>>()?;
                            for j in freq_mc.clone() {
                                if j.is_dir(){
                                    let file_mc = read_dir(j)?
                                        .map(|res| res.map(|a| a.path()))
                                        .collect::<Result<Vec<_>, io::Error>>()?;
                                    for fmc in file_mc {
                                        mc_vec.push(fmc);
                                    }
                                }
                            }
                }
            }
        }

        let data_struct = DirStruct {
            _outer_directory,
            mc_path: mc_vec,
            _mc_samples,
            _start_date: String::from(_start_date),
            _observation_time: String::from(obs_time),
            entry_point: String::from(
                            match path.to_str() {
                                Some(v) => v,
                                None => "Empty",
                            }
                        ),
            _phantom: &PhantomData::<f32>,
        };

        Ok(data_struct)
    } // end of the new function
}

impl <'a> DirStruct <'a> {

    // We have to reimplement this functions in order to load the files in parallel 
    // using the full disk bandwidth

    pub fn create_observations(&self, id_mc_ref: &str, sky_t: Vec<f32>) -> Obs {
        
        let mut tod_vec: Vec<Vec<f32>> = Vec::new();
        let mut pix_vec: Vec<Vec<i32>> = Vec::new();
        let mut det_vec: Vec<String>   = Vec::new();

        let mc_id =   match id_mc_ref.split("43").skip(1).next() {
                                            Some(c) => match c.split(".").next(){
                                                Some(m_o) => m_o,
                                                None => "Primo Err.",
                                            },
                                            None => "Secondo Err.",
                                        };
        let mc_id_u8: u8 = String::from(mc_id).parse::<u8>().unwrap();
 
        let total_file = self.get_mc_path();
        let total_file_string= total_file.iter()
                                                    .map(|f| match f.to_str() {Some(p) => p, None=> ""});
        let files = total_file_string
                                                    .filter(
                                                        |&s| wildmatch::WildMatch::new(format!("*{}*",id_mc_ref).as_str())
                                                        .matches(s)
                                                    );

        for i in files { // i -> same Monte Carlo iteration

            let mut a = NpzReader::new(File::open(i).unwrap()).unwrap();
            let c: Array1<f32> = a.by_name("arr_0.npy").unwrap();
            tod_vec.push(c.to_vec());
            /************************************************************************************ */
     

            /*************/ // Extract the detector name
            let det_i = match 
                                match i.split("/")
                                      .last() {
                                          Some(det) => det, 
                                          None=>""
                                        }.split(".").next(){
                                            Some(det) => det, 
                                            None=>""
                                        };

            det_vec.push(String::from(det_i));
            /**************/ // det_i is the name of the detector
            // Load pixel stream

            let part1 = match i.split(self.entry_point.as_str()).skip(1).next() {Some(v) => v, None => "Err."}; 
            let part2 = match part1.split("/").next() {Some(v) => v, None => "Err."}; // giorno 2022....

            let _pix_file_name = format!("{}{}/{}_pix.npz", self.entry_point, part2, det_i);
            
            let mut _pix_point: Vec<i32> =Vec::new();
            let pix_path = Path::new(_pix_file_name.as_str());

            // MEMORY ALLOCATION!!!!!!!
            let mut a = NpzReader::new(File::open(pix_path).unwrap()).unwrap();
            let c: Array1<i32> = a.by_name("arr_0.npy").unwrap();
            pix_vec.push(c.to_vec());
            

            let pix_name_string = match pix_path.to_str() {Some(p)=> p, None=> ""};
            println!("TOD: {}, PIX: {}", i.bright_green(), pix_name_string.bright_green());

        } // end of the main loop scope

        /* The signature of the new function of the Obs struct. The struc is constructed using 
           a &'static Vec<Vec<f32>> type for tods container. In this way I can pass it easy to the various threads
           that create the maps for each detector */
     
        let obs = Obs::new(
            String::from("Start"), 
            String::from("Stop"), 
            det_vec, 
            mc_id_u8 as u8, 
            1.0, 
            0.01,
            tod_vec,
            sky_t,
            pix_vec,
        );
        println!("{}", "OBS BUILT".bright_green());

        obs
    }
}

impl <'a> DirStruct <'a> {
    pub fn _get_outer_directory(&self) -> Vec<PathBuf> {
        self._outer_directory.clone()
    }
    pub fn _get_mc_samples(&self) -> usize {
        self._mc_samples
    }

    pub fn _get_obs_time(&self) -> String {
        self._observation_time.clone()
    }

    pub fn get_mc_path(&self) -> Vec<PathBuf> {
        self.mc_path.clone()
    }

}
