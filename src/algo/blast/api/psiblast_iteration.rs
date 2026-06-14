use std::collections::BTreeSet;

pub struct CPsiBlastIterationState {
    pub total_num_iterations_to_do: u32,
    pub iterations_done: u32,
    pub previous_data: BTreeSet<String>,
    pub current_data: BTreeSet<String>,
}
