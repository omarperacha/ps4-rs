use std::fs;
use std::path::Path;
use std::sync::{Mutex};
use std::collections::{HashSet, HashMap};
use std::cmp;
use csv::WriterBuilder;
use serde::Serialize;

use crate::common::levenshtein::*;

pub fn get_input_seqs(in_path: String, out_path: String) -> bool {

    let pool = rayon::ThreadPoolBuilder::new().num_threads(8).build().unwrap();

    let dir = &Path::new(&in_path);
    let res_all = Mutex::new(HashMap::new());
    let ss_all = Mutex::new(HashMap::new());
    let first_res_num_all = Mutex::new(HashMap::new());
    let entry_count = Mutex::new(0);

    pool.scope(|s| {
        for entry in fs::read_dir(dir).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path().into_os_string().into_string().unwrap();
            let file = path.split("/").last().unwrap().to_owned();
            let pdb_code = file.split(".").next().unwrap().to_owned();
            let res_all = &res_all;
            let ss_all = &ss_all;
            let first_res_num_all = &first_res_num_all;
            let entry_count = &entry_count;

            s.spawn(move | _s| {
                let data_vec_vec = get_dssp_data(&path, &pdb_code);
                for data_vec in data_vec_vec {
                    
                    // vector of strings
                    let sets_to_dedupe = vec![
                        "ps4_data/data/cb513/CB513_HHblits.csv",
                        "ps4_data/data/data.csv"
                    ];

                    let mut f = first_res_num_all.lock().unwrap();
                    let first_res = data_vec.1;

                    let mut r = res_all.lock().unwrap();
                    let res = data_vec.2;

                    let mut s = ss_all.lock().unwrap();
                    let ss = data_vec.3;

                    let mut e = entry_count.lock().unwrap();

                    let res_string: String = res.iter().collect();

                    if res.len() < 16 {
                        println!("RES TOO SHORT OR MISSING");
                        continue;
                    }

                    let mut dedupes_passed = 0;
                    let num_sets = sets_to_dedupe.len();

                    //loop through sets to dedupe
                    for set in sets_to_dedupe {
                        if is_above_thresh_with_set(&res_string, set, get_similarity_thresh()) {
                            dedupes_passed += 1;
                        } else {
                            println!("BELOW THRESH FOR {}", set);
                        }
                    }

                    if dedupes_passed == num_sets {
                        
                        *e += 1;
                        println!("{}: {} - seq len {}", *e, data_vec.0, ss.len());
                        
                        f.insert(data_vec.0.clone(), first_res);
                        r.insert(data_vec.0.clone(), res);
                        s.insert(data_vec.0, ss);

                    } 
                }
            }); 

        }
    });

    let res_unlocked = res_all.into_inner().unwrap();
    let ss_unlocked = ss_all.into_inner().unwrap();
    let first_res_unlocked = first_res_num_all.into_inner().unwrap();

    let deduped_chains: HashSet<String> = dedupe_chains(res_unlocked.clone());

    println!("\nDEDUPED CHAINS FINAL COUNT: {}", deduped_chains.len());

    let serialized = serde_pickle::to_vec(&deduped_chains, Default::default()).unwrap();
    fs::write("./res/final_deduped_filter.pkl", serialized).expect("Unable to write file");

    write_csv(&out_path, deduped_chains, res_unlocked, ss_unlocked, first_res_unlocked);

    true

}

// MARK: - RESIDUE data
fn get_dssp_data(dssp: &str, pdb_code: &str) -> Vec<(String, i32, Vec<char>, Vec<char>)> {
    let data = fs::read_to_string(dssp).expect("Unable to read file");

    let mut post = Vec::new();

    let mut res_data = Vec::new();
    let mut ss_data = Vec::new();

    let mut chain = ' ';
    let mut current_res = -999;
    let mut first_res = -999;
    let mut residues_found = false;

    for l in data.lines() {
        if residues_found {

            if l.contains("!") {
                continue;
            }

            if chain != get_chain_for(l) && chain != ' ' {
                post.push((pdb_code.to_owned() + &chain.to_string(), first_res, res_data, ss_data));
                res_data = Vec::new();
                ss_data = Vec::new();
                current_res = -999;
                first_res = -999;
            }

            chain = get_chain_for(l);

            let (num, res) = get_res_num_name(l);

            if current_res > -999 && num != current_res + 1 {
                res_data = Vec::new();
                ss_data = Vec::new();
                current_res = -998; // this will allow skipping until next chain found
                first_res = -998;
                continue;
            }

            let ss = get_ss_for(l);

            first_res = if first_res == -999 { num } else { first_res };
            
            current_res = num;
            res_data.push(res);
            ss_data.push(ss);
        } else {
            residues_found = l.split_whitespace().next().unwrap() == "#";
        }
    }
    
    post.push((pdb_code.to_owned() + &chain.to_string(), first_res, res_data, ss_data));

    return post;
}

fn get_res_num_name(residue: &str) -> (i32, char) {
    let res_id_str: String = residue[5..10].split_whitespace().collect();
    let res_name: u8 = residue.as_bytes()[13];

    (res_id_str.parse::<i32>().unwrap(), res_name as char)
}

fn get_chain_for(residue: &str) -> char {
    let chain: u8 = residue.as_bytes()[11];
    chain as char
}

fn get_ss_for(residue: &str) -> char {
    let ss: u8 = residue.as_bytes()[16];
    let ss_char = ss as char;
    if ss_char == ' ' { return 'C' } else { return ss_char };
}

fn get_beta_data_for(residue: &str) -> (i32, i32, char) {
    let res_id_str: String = residue[5..10].split_whitespace().collect();
    let res_id = res_id_str.parse::<i32>().unwrap();

    let beta_pair_str: String = residue[29..33].split_whitespace().collect();
    let beta_pair = beta_pair_str.parse::<i32>().unwrap();

    let beta_id: u8 = residue.as_bytes()[33];
    (res_id, beta_pair, beta_id as char)
}


// MARK: - Deduplication
pub fn dedupe_chains(chain_map: HashMap<String, Vec<char>>) -> HashSet<String> {

    let num_threads = 8;

    let pool = rayon::ThreadPoolBuilder::new().num_threads(num_threads).build().unwrap();
    let deduped_chunks: Mutex<Vec<HashSet<String>>> = Mutex::new(Vec::new());

    let chains: Vec<String> = chain_map.keys().cloned().collect();
    let interval = chain_map.len() / (num_threads - 1);
    let mut in_seqs = Vec::new();

    for i in 0..num_threads {
        let in_seq = &chains[i*interval..cmp::min((i+1)*interval, chains.len())];
        let in_vec = in_seq.to_vec();
        in_seqs.push(in_vec);
    }


    pool.scope(|s| {

        let deduped_chunks = &deduped_chunks;

        for i in 0..num_threads {

            let chain_map_safe = chain_map.clone();
            let in_seqs_safe = in_seqs.clone();
            s.spawn(move | _s| {
                let is = in_seqs_safe;
                let in_seqs_chunk = &*is;
                let in_keys = &in_seqs_chunk[i];

                let mut deduped: HashSet<String> = HashSet::new();
                    
                for k in in_keys {
                    let mut insert_to_deduped = true;

                    let chars = &chain_map_safe[k];
                    let in_string: String = chars.iter().collect();

                    for j in deduped.clone() {
                        let compare_string: String = chain_map_safe[&j].iter().collect();
                        if !is_above_thresh_with_seq(&in_string, &compare_string, get_similarity_thresh()) {
                            insert_to_deduped = false;
                            println!("BELOW INTRA-SET THRESH");
                            break;
                        }
                    }
        
                    if insert_to_deduped {
                        deduped.insert(k.to_string());
                        println!("DEDUPED CHAINS {} COUNT: {}", i,  deduped.len());
                    }
                }

                let mut d = deduped_chunks.lock().unwrap();
                d.push(deduped);
            });
        }
    });

    let mut deduped_chunks = deduped_chunks.into_inner().unwrap();

    while deduped_chunks.len() > 1 {

        println!("\n DIVIDE AND CONQUER, N = {}", deduped_chunks.len());

        let interval = deduped_chunks.len() / 2;

        let temp_deduped_chunks: Mutex<Vec<HashSet<String>>> = Mutex::new(Vec::new());

        pool.scope(|s| {

            let temp_deduped_chunks = &temp_deduped_chunks;
        
            for i in 0..interval {

                let chains = chain_map.clone();
                let set1 = &deduped_chunks[i * 2];
                let set2 = &deduped_chunks[(i * 2) + 1];

                s.spawn(move | _s| {
                    let filtered_chunk = dedupe_sets(set1, set2, chains);
                    let mut t = temp_deduped_chunks.lock().unwrap();
                    t.push(filtered_chunk);
                });
            }
        });

        deduped_chunks = temp_deduped_chunks.into_inner().unwrap();
    }

    let final_deduped = deduped_chunks[0].to_owned();
    final_deduped

}


fn dedupe_sets(set1: &HashSet<String>, set2: &HashSet<String>, chains: HashMap<String, Vec<char>>) -> HashSet<String> {
    
    let pool = rayon::ThreadPoolBuilder::new().num_threads(8).build().unwrap();

    let deduped: Mutex<HashSet<String>> = Mutex::new(HashSet::new());

    let chains_safe = Mutex::new(chains.clone());
    let chains_safe_2 = Mutex::new(chains.clone());

    pool.scope(|s| {

        let deduped = &deduped;
        let chains_safe = &chains_safe;
        let chains_safe_2 = &chains_safe_2;

        for k in set1 {

            let c = chains_safe.lock().unwrap();
            let in_chars = &c[k];
            let in_string: String = in_chars.to_owned().iter().collect();

            s.spawn(move | _s| {

                let mut insert_to_deduped = true;

                for j in set2 {
                    
                    let c2 = chains_safe_2.lock().unwrap();
                    let compare_chars = &c2[j];
                    let compare_string: String = compare_chars.to_owned().iter().collect();
                    if !is_above_thresh_with_seq(&in_string, &compare_string, get_similarity_thresh()) {
                        insert_to_deduped = false;
                        println!("BELOW INTER-SET THRESH");
                        break;
                    }
                }

                if insert_to_deduped {
                    let mut d = deduped.lock().unwrap();
                    d.insert(k.to_string());
                    println!("FILLTERED CHAINS COUNT: {}", d.len());
                }

            });

        }
    });

    let mut final_deduped = deduped.into_inner().unwrap();
    final_deduped.extend(set2.to_owned());

    final_deduped
}

fn get_similarity_thresh() -> Option<f32> {
    Some(60.0)
}

// MARK: - CSV writing
#[derive(Serialize)]
struct Row<'a> {
    chain_id: &'a str,
    first_res: i32,
    input: &'a str,
    dssp8: &'a str,
}

fn write_csv(
    out_path: &str,
    chains: HashSet<String>, 
    res: HashMap<String, Vec<char>>, 
    ss: HashMap<String, Vec<char>>, 
    first_res: HashMap<String, i32>) 
    {

    let path = out_path;
    if Path::new(path).exists() {
        fs::remove_file(path).unwrap();
    }

    fs::File::create(path).unwrap();

    let mut wtr = WriterBuilder::new().from_path(path).unwrap();

    for k in chains {

        let res_string: String = res[&k].iter().collect();
        let ss_string: String = ss[&k].iter().collect();

        assert_eq!(res_string.len(), ss_string.len());

        wtr.serialize(Row {
            chain_id: &k,
            first_res: first_res[&k],
            input: &res_string,
            dssp8: &ss_string,
        }).unwrap();
    }

    wtr.flush().unwrap();
}