use std::fs::File;
use std::sync::Mutex;
use edit_distance::edit_distance;

pub fn compare_sets(in_set: &str, compare_set: &str, threshold: Option<f32>) -> bool {

    let pool = rayon::ThreadPoolBuilder::new().num_threads(8).build().unwrap();

    let thresh = threshold.unwrap_or(40.0);

    let mut i_df = load_csv(in_set);

    let dists: Mutex<Vec<f32>> =  Mutex::new(Vec::new());

    pool.scope(|s| {

        let dists = &dists;
        
        
        for i_result in i_df.records() {

            s.spawn(move | _s| {
                let mut i_column = -1;
                let mut i_input_column_idx = -1;
                let i_record = i_result.unwrap();
                let mut c_df = load_csv(compare_set);

                let mut min_dist = 100.0;

                for i_str in i_record.iter() {
                    i_column += 1;

                    if i_str == "input" {
                        i_input_column_idx = i_column;
                    }

                    if i_input_column_idx != i_column {
                        continue;
                    }

                    for c_result in c_df.records() {
                        let mut c_column = -1;
                        let mut c_input_column_idx = -1;
                        let c_record = c_result.unwrap();
                        for c_str in c_record.iter() {
                            c_column += 1;

                            if c_str == "input" {
                                c_input_column_idx = c_column;
                            }

                            if c_input_column_idx != c_column {
                                continue;
                            }

                            let i_clean = i_str.replace(" ", "");
                            let c_clean = c_str.replace(" ", "");

                            let dist = edit_distance(&i_clean, &c_clean);
                            let pct = 100.0 * dist as f32 / i_clean.len() as f32;
                            min_dist = if pct < min_dist {pct} else {min_dist};
                            break;
                        }
                        if min_dist < thresh { 
                            min_dist = 0.0;
                            break 
                        };
                    }
                    let mut d = dists.lock().unwrap();
                    println!("i: {}", d.len());
                    println!("\t{}", min_dist);
                    d.push(min_dist);
                    break;
                }
            });
        }
    });

    let (mean_all, pct_all) = mean_above_thresh(dists.into_inner().unwrap());

    println!("\nDONE - pct below: {}, mean above: {}", pct_all, mean_all);

    true

}

pub fn is_above_thresh_with_set(in_seq: &str, compare_set: &str, threshold: Option<f32>) -> bool {
    
    let mut c_df = load_csv(compare_set);

    let mut min_dist = 100.0;

    let thresh = threshold.unwrap_or(20.0);

    for c_result in c_df.records() {
        let mut c_column = -1;
        let mut c_input_column_idx = -1;
        let c_record = c_result.unwrap();
        for c_str in c_record.iter() {
            
            c_column += 1;

            if c_str == "input" {
                c_input_column_idx = c_column;
            }

            if c_input_column_idx != c_column {
                continue;
            }

            let i_clean = in_seq.replace(" ", "");
            let c_clean = c_str.replace(" ", "");

            let dist = edit_distance(&i_clean, &c_clean);
            let pct = 100.0 * dist as f32 / i_clean.len() as f32;
            min_dist = if pct < min_dist {pct} else {min_dist};
            break;
        }
        if min_dist < thresh { 
            return false; 
        };
    }
    return true;
}


pub fn is_above_thresh_with_seq(in_seq: &str, compare_seq: &str, threshold: Option<f32>) -> bool {
    
    let thresh = threshold.unwrap_or(20.0);

    let i_clean = in_seq.replace(" ", "");
    let c_clean = compare_seq.replace(" ", "");

    let dist = edit_distance(&i_clean, &c_clean);
    let pct = 100.0 * dist as f32 / i_clean.len() as f32;

    return pct > thresh

}

// MARK: - Private
fn load_csv(path: &str) -> csv::Reader<File> {
    let file = File::open(format!("{}", path)).unwrap();
    let rdr = csv::Reader::from_reader(file);

    return rdr
}

fn mean_above_thresh(vec: Vec<f32>) -> (f32, f32) {
    let mut sum = 0.0;
    let mut count = 0.0;

    for i in &vec {
        if i > &0.0 {
            sum += i;
            count += 1.0;
        }
    }
    (sum / count, 100.0 * (1.0 - (count / vec.len() as f32)) )
}
