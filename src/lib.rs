use pyo3::prelude::*;
use dssp::parsing::get_input_seqs;
use common::levenshtein::compare_sets;
mod dssp;
mod common;

#[pyfunction]
fn parse_new_samples(in_path: String, out_path: String) -> PyResult<bool> {
    Ok(get_input_seqs(in_path, out_path))
}

#[pyfunction]
fn compare_two_sets(set_path_a: String, set_path_b: String) -> PyResult<bool> {
    Ok(compare_sets(&set_path_a, &set_path_b, None))
}

/// A Python module implemented in Rust.
#[pymodule]
fn ps4_rs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_new_samples, m)?)?;
    m.add_function(wrap_pyfunction!(compare_two_sets, m)?)?;
    Ok(())
}