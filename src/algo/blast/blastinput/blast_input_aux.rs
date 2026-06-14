// Upstream sources:
// - ncbi-blast-2.17.0+-src/c++/include/algo/blast/blastinput/blast_input_aux.hpp
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/blastinput/blast_input_aux.cpp

use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AutoOutputFileReset {
    pub file_name: String,
    pub file_stream: Option<String>,
    pub version: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ArgAllowMaximumFileNameLength {
    pub max_length: u32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArgAllowValuesGreaterThanOrEqual {
    pub min_value: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArgAllowValuesLessThanOrEqual {
    pub max_value: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArgAllowValuesBetween {
    pub min_value: f64,
    pub max_value: f64,
    pub inclusive: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArgAllowIntegerSet {
    pub allowed_values: BTreeSet<i64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ArgAllowDoubleSet {
    pub allowed_values: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArgAllowStringSet {
    pub allowed_values: BTreeSet<String>,
}
