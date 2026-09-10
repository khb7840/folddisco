use std::collections::HashMap;
use std::fs;
use std::process;

use folddisco::controller::io::read_structure_from_path;
use folddisco::structure::coordinate::{approx_cb, Coordinate};
use folddisco::structure::core::Structure;
use pico_args::Arguments;

#[derive(Debug, Clone)]
struct ResidueRecord {
    chain: u8,
    serial: u64,
    res_name: [u8; 3],
    n: Option<Coordinate>,
    ca: Option<Coordinate>,
    cb: Option<Coordinate>,
    c: Option<Coordinate>,
}

#[derive(Debug, Clone, PartialEq)]
struct InterfaceResidue {
    chain: u8,
    serial: u64,
    min_distance: f32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum SelectionMode {
    Top(usize),
    Threshold(f32),
}

fn main() {
    if let Err(err) = run() {
        eprintln!("{err}");
        process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let mut args = Arguments::from_env();

    if args.contains(["-h", "--help"]) {
        print_usage();
        return Ok(());
    }

    let input_path: String = args
        .value_from_str(["-p", "--path"])
        .map_err(|_| "Missing required argument: --path <PDB|mmCIF path>".to_string())?;
    let chain_pair_raw: String = args
        .value_from_str(["-c", "--chains"])
        .map_err(|_| "Missing required argument: --chains <CHAIN1,CHAIN2>".to_string())?;
    let output_path: Option<String> = args.opt_value_from_str(["-o", "--output"]).map_err(|e| e.to_string())?;
    let top_n: Option<usize> = args.opt_value_from_str(["-n", "--top"]).map_err(|e| e.to_string())?;
    let threshold: Option<f32> = args
        .opt_value_from_str(["-d", "--distance-threshold"])
        .map_err(|e| e.to_string())?;

    let remaining = args.finish();
    if !remaining.is_empty() {
        return Err(format!("Unexpected arguments: {:?}", remaining));
    }

    let selection_mode = match (top_n, threshold) {
        (Some(_), Some(_)) => {
            return Err("Use either --top or --distance-threshold, not both".to_string())
        }
        (Some(0), None) => return Err("--top must be greater than zero".to_string()),
        (Some(n), None) => SelectionMode::Top(n),
        (None, Some(value)) if value < 0.0 => {
            return Err("--distance-threshold must be non-negative".to_string())
        }
        (None, Some(value)) => SelectionMode::Threshold(value),
        (None, None) => {
            return Err("One of --top <N> or --distance-threshold <FLOAT> is required".to_string())
        }
    };

    let (chain_a, chain_b) = parse_chain_pair(&chain_pair_raw)?;
    let structure = read_structure_from_path(&input_path)
        .ok_or_else(|| format!("Failed to read structure from {}", input_path))?;
    let residues = collect_residues(&structure);
    let selected = extract_interface_residues(&residues, chain_a, chain_b, selection_mode)?;
    let query_line = format_query_line(&input_path, &selected);

    if let Some(path) = output_path {
        fs::write(&path, format!("{query_line}\n"))
            .map_err(|e| format!("Failed to write output to {path}: {e}"))?;
    } else {
        println!("{query_line}");
    }

    Ok(())
}

fn print_usage() {
    eprintln!(
        "Extract inter-chain interface residues and print them in Folddisco query-file form.\n\
         \n\
         Usage:\n\
           cargo run --bin extract_interface_query -- --path <STRUCTURE> --chains <A,B> --top <N>\n\
           cargo run --bin extract_interface_query -- --path <STRUCTURE> --chains <A,B> --distance-threshold <ANGSTROM>\n\
         \n\
         Options:\n\
           -p, --path <PATH>                  Input PDB/mmCIF file\n\
           -c, --chains <CHAIN1,CHAIN2>       Chain pair to analyze\n\
           -n, --top <INT>                    Keep the top N interface residues by minimum inter-chain distance\n\
           -d, --distance-threshold <FLOAT>   Keep residues whose minimum inter-chain distance is within the threshold\n\
           -o, --output <PATH>                Write the Folddisco query line to a file instead of stdout\n\
           -h, --help                         Show this help message"
    );
}

fn parse_chain_pair(raw: &str) -> Result<(u8, u8), String> {
    let mut parts = raw.split(',').map(str::trim);
    let left = parts
        .next()
        .ok_or_else(|| "Missing first chain in --chains".to_string())?;
    let right = parts
        .next()
        .ok_or_else(|| "Missing second chain in --chains".to_string())?;

    if parts.next().is_some() {
        return Err("--chains must contain exactly two comma-separated chain identifiers".to_string());
    }

    let left = parse_single_chain(left)?;
    let right = parse_single_chain(right)?;

    if left == right {
        return Err("--chains requires two different chains".to_string());
    }

    Ok((left, right))
}

fn parse_single_chain(raw: &str) -> Result<u8, String> {
    if raw.len() != 1 {
        return Err(format!(
            "Chain identifier '{raw}' is invalid; only single-character chain IDs are supported"
        ));
    }
    Ok(raw.as_bytes()[0])
}

fn collect_residues(structure: &Structure) -> Vec<ResidueRecord> {
    let atoms = &structure.atom_vector;
    let mut residues = Vec::new();
    let mut current: Option<ResidueRecord> = None;

    for idx in 0..atoms.len() {
        let atom = atoms.get(idx);
        let changed = match &current {
            Some(record) => record.chain != atom.chain || record.serial != atom.res_serial,
            None => true,
        };

        if changed {
            if let Some(record) = current.take() {
                residues.push(record);
            }
            current = Some(ResidueRecord {
                chain: atom.chain,
                serial: atom.res_serial,
                res_name: atom.res_name,
                n: None,
                ca: None,
                cb: None,
                c: None,
            });
        }

        if let Some(record) = current.as_mut() {
            let coord = atom.get_coordinate();
            match atom.atom_name.as_slice() {
                b" N  " => record.n = Some(coord),
                b" CA " => record.ca = Some(coord),
                b" CB " => record.cb = Some(coord),
                b" C  " => record.c = Some(coord),
                _ => {}
            }
        }
    }

    if let Some(record) = current {
        residues.push(record);
    }

    residues
}

fn extract_interface_residues(
    residues: &[ResidueRecord],
    chain_a: u8,
    chain_b: u8,
    selection_mode: SelectionMode,
) -> Result<Vec<InterfaceResidue>, String> {
    let left: Vec<&ResidueRecord> = residues.iter().filter(|r| r.chain == chain_a).collect();
    let right: Vec<&ResidueRecord> = residues.iter().filter(|r| r.chain == chain_b).collect();

    if left.is_empty() {
        return Err(format!("No residues found for chain {}", chain_a as char));
    }
    if right.is_empty() {
        return Err(format!("No residues found for chain {}", chain_b as char));
    }

    let left_best = min_partner_distance_map(&left, &right);
    let right_best = min_partner_distance_map(&right, &left);

    let mut combined = Vec::with_capacity(left_best.len() + right_best.len());
    combined.extend(map_to_interface_residues(chain_a, &left_best));
    combined.extend(map_to_interface_residues(chain_b, &right_best));
    combined.sort_by(|a, b| {
        a.min_distance
            .total_cmp(&b.min_distance)
            .then_with(|| a.chain.cmp(&b.chain))
            .then_with(|| a.serial.cmp(&b.serial))
    });

    let selected = match selection_mode {
        SelectionMode::Top(n) => combined.into_iter().take(n).collect(),
        SelectionMode::Threshold(threshold) => combined
            .into_iter()
            .filter(|residue| residue.min_distance <= threshold)
            .collect(),
    };

    Ok(selected)
}

fn min_partner_distance_map(
    source: &[&ResidueRecord],
    partner: &[&ResidueRecord],
) -> HashMap<u64, f32> {
    let mut best = HashMap::new();

    for residue in source {
        let Some(residue_coord) = residue.reference_coordinate() else {
            continue;
        };
        let mut min_distance: Option<f32> = None;

        for partner_residue in partner {
            let Some(partner_coord) = partner_residue.reference_coordinate() else {
                continue;
            };
            let distance = residue_coord.calc_distance(&partner_coord);
            min_distance = Some(match min_distance {
                Some(current) => current.min(distance),
                None => distance,
            });
        }

        if let Some(distance) = min_distance {
            best.insert(residue.serial, distance);
        }
    }

    best
}

fn map_to_interface_residues(chain: u8, distances: &HashMap<u64, f32>) -> Vec<InterfaceResidue> {
    distances
        .iter()
        .map(|(serial, min_distance)| InterfaceResidue {
            chain,
            serial: *serial,
            min_distance: *min_distance,
        })
        .collect()
}

fn format_query_line(input_path: &str, residues: &[InterfaceResidue]) -> String {
    let residues = residues
        .iter()
        .map(|residue| format!("{}{}", residue.chain as char, residue.serial))
        .collect::<Vec<_>>()
        .join(",");
    format!("{input_path}\t{residues}")
}

impl ResidueRecord {
    fn reference_coordinate(&self) -> Option<Coordinate> {
        if let Some(cb) = self.cb {
            return Some(cb);
        }
        if self.res_name == *b"GLY" {
            if let (Some(ca), Some(n), Some(c)) = (self.ca, self.n, self.c) {
                return Some(approx_cb(&ca, &n, &c));
            }
        }
        self.ca
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn residue(chain: u8, serial: u64, res_name: [u8; 3], ca: [f32; 3], cb: Option<[f32; 3]>) -> ResidueRecord {
        ResidueRecord {
            chain,
            serial,
            res_name,
            n: Some(Coordinate::new(ca[0] - 1.0, ca[1], ca[2])),
            ca: Some(Coordinate::new(ca[0], ca[1], ca[2])),
            cb: cb.map(|coord| Coordinate::new(coord[0], coord[1], coord[2])),
            c: Some(Coordinate::new(ca[0] + 1.0, ca[1], ca[2])),
        }
    }

    #[test]
    fn parses_chain_pair() {
        assert_eq!(parse_chain_pair("A,B").unwrap(), (b'A', b'B'));
        assert!(parse_chain_pair("A").is_err());
        assert!(parse_chain_pair("AA,B").is_err());
        assert!(parse_chain_pair("A,A").is_err());
    }

    #[test]
    fn sorts_top_residues_by_minimum_distance() {
        let residues = vec![
            residue(b'A', 1, *b"ALA", [0.0, 0.0, 0.0], Some([0.0, 0.0, 0.0])),
            residue(b'A', 2, *b"ALA", [5.0, 0.0, 0.0], Some([5.0, 0.0, 0.0])),
            residue(b'B', 10, *b"ALA", [1.0, 0.0, 0.0], Some([1.0, 0.0, 0.0])),
            residue(b'B', 11, *b"ALA", [9.0, 0.0, 0.0], Some([9.0, 0.0, 0.0])),
        ];

        let selected = extract_interface_residues(&residues, b'A', b'B', SelectionMode::Top(3)).unwrap();

        assert_eq!(selected.len(), 3);
        assert!(selected[0].min_distance <= selected[1].min_distance);
        assert!(selected[1].min_distance <= selected[2].min_distance);
        assert_eq!((selected[0].chain, selected[0].serial), (b'A', 1));
        assert_eq!((selected[1].chain, selected[1].serial), (b'B', 10));
    }

    #[test]
    fn filters_threshold_residues() {
        let residues = vec![
            residue(b'A', 1, *b"ALA", [0.0, 0.0, 0.0], Some([0.0, 0.0, 0.0])),
            residue(b'A', 2, *b"ALA", [7.0, 0.0, 0.0], Some([7.0, 0.0, 0.0])),
            residue(b'B', 10, *b"ALA", [2.0, 0.0, 0.0], Some([2.0, 0.0, 0.0])),
            residue(b'B', 11, *b"ALA", [20.0, 0.0, 0.0], Some([20.0, 0.0, 0.0])),
        ];

        let selected =
            extract_interface_residues(&residues, b'A', b'B', SelectionMode::Threshold(3.0)).unwrap();

        assert_eq!(
            selected
                .iter()
                .map(|residue| (residue.chain, residue.serial))
                .collect::<Vec<_>>(),
            vec![(b'A', 1), (b'B', 10)]
        );
    }

    #[test]
    fn formats_query_line() {
        let formatted = format_query_line(
            "/tmp/example.pdb",
            &[
                InterfaceResidue {
                    chain: b'A',
                    serial: 1,
                    min_distance: 1.0,
                },
                InterfaceResidue {
                    chain: b'B',
                    serial: 5,
                    min_distance: 2.0,
                },
            ],
        );
        assert_eq!(formatted, "/tmp/example.pdb\tA1,B5");
    }

    #[test]
    fn glycine_without_cb_uses_approximation() {
        let gly = residue(b'A', 1, *b"GLY", [0.0, 0.0, 0.0], None);
        assert!(gly.reference_coordinate().is_some());
    }

    #[test]
    fn real_structure_returns_query_line() {
        let input_path = "/home/runner/work/folddisco/folddisco/query/4CHA.pdb";
        let structure = read_structure_from_path(input_path).unwrap();
        let residues = collect_residues(&structure);
        let selected =
            extract_interface_residues(&residues, b'B', b'C', SelectionMode::Top(6)).unwrap();
        let line = format_query_line(input_path, &selected);

        assert!(line.starts_with("/home/runner/work/folddisco/folddisco/query/4CHA.pdb\t"));
        assert!(line.split('\t').nth(1).is_some_and(|value| !value.is_empty()));
        assert!(selected.len() <= 6);
    }
}
