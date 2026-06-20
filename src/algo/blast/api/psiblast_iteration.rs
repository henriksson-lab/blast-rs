use std::collections::BTreeSet;

pub struct CPsiBlastIterationState {
    pub total_num_iterations_to_do: u32,
    pub iterations_done: u32,
    pub previous_data: BTreeSet<String>,
    pub current_data: BTreeSet<String>,
}

impl CPsiBlastIterationState {
    pub fn new(num_iterations: u32) -> Self {
        Self {
            total_num_iterations_to_do: num_iterations,
            iterations_done: 0,
            previous_data: BTreeSet::new(),
            current_data: BTreeSet::new(),
        }
    }

    pub fn has_more_iterations(&self) -> bool {
        if self.total_num_iterations_to_do == 0
            || self.iterations_done < self.total_num_iterations_to_do
        {
            return true;
        }
        false
    }

    pub fn has_converged(&self) -> bool {
        if self.iterations_done <= 1 {
            return false;
        }

        if !self.previous_data.is_empty() && self.current_data.is_empty() {
            return true;
        }

        let mut retval = true;
        for current_id in &self.current_data {
            if !self.previous_data.contains(current_id) {
                retval = false;
                break;
            }
        }
        retval
    }

    pub fn is_active(&self) -> bool {
        self.has_more_iterations() && !self.has_converged()
    }

    pub fn x_throw_exception_on_logic_error(&self) {
        if !self.is_active() {
            panic!(
                "Should not modify a PSI-BLAST iteration after it has converged or exhausted its iterations"
            );
        }
    }

    pub fn get_previously_found_seq_ids(&self) -> BTreeSet<String> {
        self.current_data.clone()
    }

    pub fn advance(&mut self, list: BTreeSet<String>) {
        self.x_throw_exception_on_logic_error();
        self.previous_data = std::mem::take(&mut self.current_data);
        self.current_data = list;
        self.iterations_done += 1;
    }

    pub fn get_iteration_number(&self) -> u32 {
        self.iterations_done + 1
    }
}
