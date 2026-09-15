//! Chain identifiers wide enough for multi-character mmCIF asym IDs.
//!
//! Chain IDs used to be a single `u8`, which silently dropped everything past
//! the first character of a multi-character mmCIF `auth_asym_id` (large cryo-EM
//! entries such as PDB 9A1O use `"10"`, `"AA"`, ...). [`ChainId`] stores up to
//! [`CHAIN_ID_MAX_LEN`] printable ASCII bytes inline, null-padded on the right,
//! so it stays `Copy` and allocation-free in the retrieval hot path where one
//! chain ID is carried per residue of every candidate structure.
//!
//! Note that [`crate::structure::atom::Atom`] deliberately keeps its
//! single-byte `chain` field: it is transmuted from Foldcomp's C `atom_t`
//! (see `Atom::from_c`), whose layout must not change. Widening happens one
//! level up, in `AtomVector` / `Structure` / `CompactStructure`, which is where
//! the CIF parser can supply the full identifier.

use std::fmt;

/// Longest chain ID that is stored in full. mmCIF asym IDs in the PDB archive
/// are at most four characters today; the extra headroom costs nothing because
/// `Option<(ChainId, u64)>` is 24 bytes either way.
pub const CHAIN_ID_MAX_LEN: usize = 8;

/// Byte substituted for anything that is not printable ASCII, so that
/// [`ChainId::as_str`] can never fail on a malformed input file.
const REPLACEMENT: u8 = b'?';

/// A chain identifier of up to [`CHAIN_ID_MAX_LEN`] printable ASCII bytes,
/// null-padded on the right.
#[derive(Copy, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct ChainId {
    bytes: [u8; CHAIN_ID_MAX_LEN],
}

impl ChainId {
    /// The empty chain ID. No parsed structure produces one, so it doubles as a
    /// "nothing seen yet" sentinel.
    pub const fn empty() -> Self {
        ChainId { bytes: [0; CHAIN_ID_MAX_LEN] }
    }

    /// Wrap a single byte, as read from column 22 of a PDB `ATOM` record.
    pub fn from_byte(byte: u8) -> Self {
        let mut bytes = [0u8; CHAIN_ID_MAX_LEN];
        bytes[0] = sanitize(byte);
        ChainId { bytes }
    }

    /// Take up to [`CHAIN_ID_MAX_LEN`] bytes. Anything longer is truncated;
    /// callers that care can compare `len()` against the input length.
    pub fn from_bytes(bytes: &[u8]) -> Self {
        let mut out = [0u8; CHAIN_ID_MAX_LEN];
        for (slot, &byte) in out.iter_mut().zip(bytes.iter()) {
            *slot = sanitize(byte);
        }
        ChainId { bytes: out }
    }

    /// Take up to [`CHAIN_ID_MAX_LEN`] bytes of a string. Anything longer is
    /// truncated.
    pub fn from_str(text: &str) -> Self {
        Self::from_bytes(text.as_bytes())
    }

    /// Number of bytes actually used, 0 for [`ChainId::empty`].
    pub fn len(&self) -> usize {
        self.bytes.iter().position(|&b| b == 0).unwrap_or(CHAIN_ID_MAX_LEN)
    }

    pub fn is_empty(&self) -> bool {
        self.bytes[0] == 0
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes[..self.len()]
    }

    pub fn as_str(&self) -> &str {
        // Every byte is forced into printable ASCII by `sanitize`, so this holds.
        std::str::from_utf8(self.as_bytes()).unwrap_or("?")
    }

    /// First byte, for the single-byte `Atom::chain` field that Foldcomp's
    /// `atom_t` layout pins down. Lossy for multi-character IDs by definition.
    pub fn first_byte(&self) -> u8 {
        self.bytes[0]
    }

    /// Whether `format!("{}{}", chain, residue)` would be ambiguous to read
    /// back, i.e. whether this chain ID needs an explicit `_` before the
    /// residue index.
    ///
    /// `A` + `21` -> `A21` parses back, because a leading letter can only be
    /// the chain. `10` + `21` -> `1021` and `AA` + `250` -> `AA250` cannot.
    pub fn needs_separator(&self) -> bool {
        self.len() != 1 || !self.bytes[0].is_ascii_alphabetic()
    }
}

/// Map anything that is not printable ASCII onto [`REPLACEMENT`], so a
/// corrupted input file cannot produce a `ChainId` that is not valid UTF-8.
fn sanitize(byte: u8) -> u8 {
    if byte.is_ascii_graphic() || byte == b' ' { byte } else { REPLACEMENT }
}

impl fmt::Display for ChainId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl fmt::Debug for ChainId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ChainId({:?})", self.as_str())
    }
}

impl From<u8> for ChainId {
    fn from(byte: u8) -> Self {
        ChainId::from_byte(byte)
    }
}

impl From<&str> for ChainId {
    fn from(text: &str) -> Self {
        ChainId::from_str(text)
    }
}

/// Shorthand for `ChainId::from_str`, handy at call sites that build chain IDs
/// from string literals, e.g. in tests.
pub fn chain(text: &str) -> ChainId {
    ChainId::from_str(text)
}

/// Whether a whole residue list has to be printed with `_` separators.
///
/// The decision is taken per list rather than per residue so that one field is
/// never half in one format and half in the other: a reader can split a field
/// on `,` and apply a single rule to every element.
pub fn residue_list_needs_separator<'a, I>(chains: I) -> bool
where
    I: IntoIterator<Item = &'a ChainId>,
{
    chains.into_iter().any(|chain| chain.needs_separator())
}

/// `A21` or `A_21`, depending on `separator`.
pub fn format_chain_residue(chain: &ChainId, residue: u64, separator: bool) -> String {
    if separator {
        format!("{}_{}", chain.as_str(), residue)
    } else {
        format!("{}{}", chain.as_str(), residue)
    }
}

/// Split one query/output token into its chain and residue parts.
///
/// Accepts, in order of precedence:
///   `AA_250` / `A_250` / `10_250` — explicit separator, any chain ID
///   `A250`                        — legacy single leading letter
///   `250`                         — no chain, caller's default applies
///
/// Returns `(None, rest)` when the token carries no chain of its own.
pub fn split_chain_and_rest(token: &str) -> (Option<ChainId>, &str) {
    if let Some(sep) = token.find('_') {
        let (chain, rest) = (&token[..sep], &token[sep + 1..]);
        // A leading `_` is the "residue not matched" placeholder in output, and
        // carries no chain; treat it like a bare residue index.
        if chain.is_empty() {
            return (None, rest);
        }
        return (Some(ChainId::from_str(chain)), rest);
    }
    match token.as_bytes().first() {
        Some(first) if first.is_ascii_alphabetic() => {
            (Some(ChainId::from_byte(*first)), &token[1..])
        }
        _ => (None, token),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_str_roundtrips() {
        for text in &["A", "AA", "10", "1", "ABCD", "AAAAAAAA", "a"] {
            assert_eq!(ChainId::from_str(text).as_str(), *text);
            assert_eq!(ChainId::from_str(text).len(), text.len());
        }
    }

    #[test]
    fn from_byte_matches_from_str() {
        assert_eq!(ChainId::from_byte(b'A'), ChainId::from_str("A"));
        assert_eq!(ChainId::from_byte(b'A').first_byte(), b'A');
    }

    #[test]
    fn empty_is_distinct_from_every_real_chain() {
        assert!(ChainId::empty().is_empty());
        assert_eq!(ChainId::empty().len(), 0);
        assert_eq!(ChainId::empty().as_str(), "");
        assert_ne!(ChainId::empty(), ChainId::from_str("A"));
        assert_ne!(ChainId::empty(), ChainId::from_byte(b' '));
    }

    #[test]
    fn over_long_ids_truncate_rather_than_panic() {
        let long = ChainId::from_str("ABCDEFGHIJ");
        assert_eq!(long.len(), CHAIN_ID_MAX_LEN);
        assert_eq!(long.as_str(), "ABCDEFGH");
    }

    #[test]
    fn non_ascii_bytes_stay_printable() {
        assert_eq!(ChainId::from_bytes(&[0xff]).as_str(), "?");
        assert_eq!(ChainId::from_bytes(&[b'A', 0x01]).as_str(), "A?");
    }

    #[test]
    fn separator_is_needed_exactly_when_concatenation_is_ambiguous() {
        assert!(!ChainId::from_str("A").needs_separator());
        assert!(!ChainId::from_str("z").needs_separator());
        assert!(ChainId::from_str("AA").needs_separator());
        assert!(ChainId::from_str("1").needs_separator());
        assert!(ChainId::from_str("10").needs_separator());
        assert!(ChainId::empty().needs_separator());
    }

    #[test]
    fn residue_list_decision_is_all_or_nothing() {
        let legacy = vec![ChainId::from_str("A"), ChainId::from_str("B")];
        assert!(!residue_list_needs_separator(legacy.iter()));
        let mixed = vec![ChainId::from_str("A"), ChainId::from_str("10")];
        assert!(residue_list_needs_separator(mixed.iter()));
    }

    #[test]
    fn formatting_follows_the_separator_flag() {
        let chain = ChainId::from_str("A");
        assert_eq!(format_chain_residue(&chain, 21, false), "A21");
        assert_eq!(format_chain_residue(&chain, 21, true), "A_21");
        assert_eq!(format_chain_residue(&ChainId::from_str("10"), 21, true), "10_21");
    }

    #[test]
    fn split_accepts_both_grammars() {
        assert_eq!(split_chain_and_rest("A250"), (Some(ChainId::from_str("A")), "250"));
        assert_eq!(split_chain_and_rest("A_250"), (Some(ChainId::from_str("A")), "250"));
        assert_eq!(split_chain_and_rest("AA_250"), (Some(ChainId::from_str("AA")), "250"));
        assert_eq!(split_chain_and_rest("10_250"), (Some(ChainId::from_str("10")), "250"));
        assert_eq!(split_chain_and_rest("250"), (None, "250"));
        assert_eq!(split_chain_and_rest("_250"), (None, "250"));
        assert_eq!(split_chain_and_rest("A250-252"), (Some(ChainId::from_str("A")), "250-252"));
        assert_eq!(split_chain_and_rest("AA_250-252"), (Some(ChainId::from_str("AA")), "250-252"));
    }
}
