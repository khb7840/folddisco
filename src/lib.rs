//! # Folddisco
//!
//! Folddisco is a tool for finding discontinuous motifs in protein structures.

pub mod cli;
pub mod controller;
pub mod geometry;
pub mod index;
pub mod structure;
pub mod utils;

/* re-export: pub use */
pub use structure::io::cif::Reader as CIFReader;
pub use structure::io::pdb::Reader as PDBReader;

pub mod prelude {
    pub use crate::measure_time;
    pub use crate::PDBReader;

    pub use crate::controller::io::{read_offset_map, save_offset_map, write_usize_vector};
    pub use crate::controller::query::{make_query_map, parse_query_string};
    pub use crate::controller::Folddisco;

    pub use crate::geometry::core::{GeometricHash, HashType};

    pub use crate::index::indextable::{load_folddisco_index, FolddiscoIndex};
    pub use crate::index::lookup::{load_lookup_from_file, save_lookup_to_file};

    pub use crate::utils::loader::load_path;
    pub use crate::utils::log::{log_msg, print_log_msg, DONE, FAIL, INFO, WARN};
    pub use crate::utils::traits::HashableSync;
}
