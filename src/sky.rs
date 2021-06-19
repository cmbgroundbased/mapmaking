use npy::NpyData;
use std::{fs::File, io::Read};

pub struct Sky {
    _nside: u32,
    sky_t: Vec<f32>,
    sky_q: Vec<f32>,
    sky_u: Vec<f32>,
}


impl Sky {
    pub fn new() -> Self {
        let mut _sky_map_t: Vec<f32> = Vec::new();
        let mut _sky_map_q: Vec<f32> = Vec::new();
        let mut _sky_map_u: Vec<f32> = Vec::new();

        let mut buffer_t = Vec::new();
        let mut buffer_q = Vec::new();
        let mut buffer_u = Vec::new();

        File::open("Sky_Maps/T_map.npy").unwrap().read_to_end(&mut buffer_t).unwrap();
        File::open("Sky_Maps/Q_map.npy").unwrap().read_to_end(&mut buffer_q).unwrap();
        File::open("Sky_Maps/U_map.npy").unwrap().read_to_end(&mut buffer_u).unwrap();
        
        _sky_map_t = NpyData::from_bytes(&mut buffer_t).unwrap().to_vec();
        _sky_map_q = NpyData::from_bytes(&mut buffer_q).unwrap().to_vec();
        _sky_map_u = NpyData::from_bytes(&mut buffer_u).unwrap().to_vec();

        Sky {
            _nside: 128,
            sky_t: _sky_map_t,
            sky_q: _sky_map_q,
            sky_u: _sky_map_u,
        }
    }
} 

impl Sky {
    pub fn get_t_sky(&self) -> Vec<f32> {
        self.sky_t.clone()
    }

    pub fn get_q_sky(&self) -> Vec<f32> {
        self.sky_q.clone()
    }

    pub fn get_u_sky(&self) -> Vec<f32> {
        self.sky_u.clone()
    }
}




