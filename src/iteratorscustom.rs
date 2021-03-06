//! Create a float iterator given the `start`, `stop` and the number of elements `steps`
pub struct FloatIterator {
    current: u32,
    current_back: u32,
    steps: u32,
    start: f32,
    end: f32,
}

impl FloatIterator {
    pub fn new(start: f32, end: f32, steps: u32) -> Self {
        FloatIterator {
            current: 0,
            current_back: steps,
            steps,
            start,
            end,
        }
    }

    pub fn length(&self) -> u32 {
        self.current - self.current_back
    }

    pub fn at(&self, pos: u32) -> f32 {
        let f_pos = pos as f32 / self.steps as f32;
        (1. - f_pos) * self.start + f_pos * self.end
    }

    fn usize_len(&self) -> usize {
        let l = self.length();
        debug_assert!(l <= ::std::usize::MAX as u32);
        l as usize
    }
}

impl Iterator for FloatIterator {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current >= self.current_back {
            return None;
        }
        let result = self.at(self.current);
        self.current += 1;
        Some(result)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let l = self.usize_len();
        (l, Some(l))
    }

    fn count(self) -> usize {
        self.usize_len()
    }
}

impl DoubleEndedIterator for FloatIterator {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.current >= self.current_back {
            return None;
        }
        self.current_back -= 1;
        let result = self.at(self.current_back);
        Some(result)
    }
}

impl ExactSizeIterator for FloatIterator {
    fn len(&self) -> usize {
        self.usize_len()
    }
}
