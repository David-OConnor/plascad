#![allow(non_snake_case)]

//! Primer melting temperature calculations. Modified from [BioPython's module here](https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/MeltingTemp.py)
//! We use an approach that calculates enthalpy and entropy of neighbors basd on empirical data,
//! and apply salt corrections based on user input concentrations of ions and primers.
//!
//! The calculations are based primarily on [SantaLucia & Hicks (2004)](https://pubmed.ncbi.nlm.nih.gov/15139820/)
//!
//! [This calculator from NorthWestern](http://biotools.nubic.northwestern.edu/OligoCalc.html) may be used
//! for QC TM, weight, and other properties. It includes detailed sources and methods.

use na_seq::{
    calc_gc,
    Nucleotide::{self, A, C, G, T},
};

use crate::primer::{IonConcentrations, MIN_PRIMER_LEN};

const R: f32 = 1.987; // Universal gas constant (Cal/C * Mol)

/// Enthalpy (dH) and entropy (dS) tables based on terminal missmatch
fn _dH_dS_tmm(nts: (Nucleotide, Nucleotide)) -> Option<(f32, f32)> {
    match nts {
        (A, A) => Some((-7.6, -21.3)),
        (A, T) => Some((-7.2, -20.4)),
        (T, A) => Some((-7.2, -21.3)),
        (C, A) => Some((-8.5, -22.7)),
        (C, G) => Some((-10.6, -27.2)),
        (G, A) => Some((-8.2, -22.2)),
        (G, C) => Some((-9.8, -24.4)),
        (G, T) => Some((-8.4, -22.4)),
        (G, G) => Some((-8.0, -19.9)),
        _ => None,
    }

    // # Terminal mismatch table (DNA)
    // # SantaLucia & Peyret (2001) Patent Application WO 01/94611
    // DNA_TMM1 = {
    //     "AA/TA": (-3.1, -7.8),
    //     "TA/AA": (-2.5, -6.3),
    //     "CA/GA": (-4.3, -10.7),
    //     "GA/CA": (-8.0, -22.5),
    //     "AC/TC": (-0.1, 0.5),
    //     "TC/AC": (-0.7, -1.3),
    //     "CC/GC": (-2.1, -5.1),
    //     "GC/CC": (-3.9, -10.6),
    //     "AG/TG": (-1.1, -2.1),
    //     "TG/AG": (-1.1, -2.7),
    //     "CG/GG": (-3.8, -9.5),
    //     "GG/CG": (-0.7, -19.2),
    //     "AT/TT": (-2.4, -6.5),
    //     "TT/AT": (-3.2, -8.9),
    //     "CT/GT": (-6.1, -16.9),
    //     "GT/CT": (-7.4, -21.2),
    //     "AA/TC": (-1.6, -4.0),
    //     "AC/TA": (-1.8, -3.8),
    //     "CA/GC": (-2.6, -5.9),
    //     "CC/GA": (-2.7, -6.0),
    //     "GA/CC": (-5.0, -13.8),
    //     "GC/CA": (-3.2, -7.1),
    //     "TA/AC": (-2.3, -5.9), "TC/AA": (-2.7, -7.0),
    //     "AC/TT": (-0.9, -1.7), "AT/TC": (-2.3, -6.3), "CC/GT": (-3.2, -8.0),
    //     "CT/GC": (-3.9, -10.6), "GC/CT": (-4.9, -13.5), "GT/CC": (-3.0, -7.8),
    //     "TC/AT": (-2.5, -6.3), "TT/AC": (-0.7, -1.2),
    //     "AA/TG": (-1.9, -4.4), "AG/TA": (-2.5, -5.9), "CA/GG": (-3.9, -9.6),
    //     "CG/GA": (-6.0, -15.5), "GA/CG": (-4.3, -11.1), "GG/CA": (-4.6, -11.4),
    //     "TA/AG": (-2.0, -4.7), "TG/AA": (-2.4, -5.8),
    //     "AG/TT": (-3.2, -8.7), "AT/TG": (-3.5, -9.4), "CG/GT": (-3.8, -9.0),
    //     "CT/GG": (-6.6, -18.7), "GG/CT": (-5.7, -15.9), "GT/CG": (-5.9, -16.1),
    //     "TG/AT": (-3.9, -10.5), "TT/AG": (-3.6, -9.8)}
}

/// Enthalpy (dH) and entropy (dS) tables based on internal missmatch
fn _dH_dS_imm(nts: (Nucleotide, Nucleotide)) -> Option<(f32, f32)> {
    match nts {
        (A, A) => Some((-7.6, -21.3)),
        (A, T) => Some((-7.2, -20.4)),
        (T, A) => Some((-7.2, -21.3)),
        (C, A) => Some((-8.5, -22.7)),
        (C, G) => Some((-10.6, -27.2)),
        (G, A) => Some((-8.2, -22.2)),
        (G, C) => Some((-9.8, -24.4)),
        (G, T) => Some((-8.4, -22.4)),
        (G, G) => Some((-8.0, -19.9)),
        _ => None,
    }

    // # Internal mismatch and inosine table (DNA)
    // # Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594
    // # Allawi & SantaLucia (1998), Biochemistry 37: 9435-9444
    // # Allawi & SantaLucia (1998), Biochemistry 37: 2170-2179
    // # Allawi & SantaLucia (1998), Nucl Acids Res 26: 2694-2701
    // # Peyret et al. (1999), Biochemistry 38: 3468-3477
    // # Watkins & SantaLucia (2005), Nucl Acids Res 33: 6258-6267
    // DNA_IMM1 = {
    //     "AG/TT": (1.0, 0.9), "AT/TG": (-2.5, -8.3), "CG/GT": (-4.1, -11.7),
    //     "CT/GG": (-2.8, -8.0), "GG/CT": (3.3, 10.4), "GG/TT": (5.8, 16.3),
    //     "GT/CG": (-4.4, -12.3), "GT/TG": (4.1, 9.5), "TG/AT": (-0.1, -1.7),
    //     "TG/GT": (-1.4, -6.2), "TT/AG": (-1.3, -5.3), "AA/TG": (-0.6, -2.3),
    //     "AG/TA": (-0.7, -2.3), "CA/GG": (-0.7, -2.3), "CG/GA": (-4.0, -13.2),
    //     "GA/CG": (-0.6, -1.0), "GG/CA": (0.5, 3.2), "TA/AG": (0.7, 0.7),
    //     "TG/AA": (3.0, 7.4),
    //     "AC/TT": (0.7, 0.2), "AT/TC": (-1.2, -6.2), "CC/GT": (-0.8, -4.5),
    //     "CT/GC": (-1.5, -6.1), "GC/CT": (2.3, 5.4), "GT/CC": (5.2, 13.5),
    //     "TC/AT": (1.2, 0.7), "TT/AC": (1.0, 0.7),
    //     "AA/TC": (2.3, 4.6), "AC/TA": (5.3, 14.6), "CA/GC": (1.9, 3.7),
    //     "CC/GA": (0.6, -0.6), "GA/CC": (5.2, 14.2), "GC/CA": (-0.7, -3.8),
    //     "TA/AC": (3.4, 8.0), "TC/AA": (7.6, 20.2),
    //     "AA/TA": (1.2, 1.7), "CA/GA": (-0.9, -4.2), "GA/CA": (-2.9, -9.8),
    //     "TA/AA": (4.7, 12.9), "AC/TC": (0.0, -4.4), "CC/GC": (-1.5, -7.2),
    //     "GC/CC": (3.6, 8.9), "TC/AC": (6.1, 16.4), "AG/TG": (-3.1, -9.5),
    //     "CG/GG": (-4.9, -15.3), "GG/CG": (-6.0, -15.8), "TG/AG": (1.6, 3.6),
    //     "AT/TT": (-2.7, -10.8), "CT/GT": (-5.0, -15.8), "GT/CT": (-2.2, -8.4),
    //     "TT/AT": (0.2, -1.5),
    //     "AI/TC": (-8.9, -25.5), "TI/AC": (-5.9, -17.4), "AC/TI": (-8.8, -25.4),
    //     "TC/AI": (-4.9, -13.9), "CI/GC": (-5.4, -13.7), "GI/CC": (-6.8, -19.1),
    //     "CC/GI": (-8.3, -23.8), "GC/CI": (-5.0, -12.6),
    //     "AI/TA": (-8.3, -25.0), "TI/AA": (-3.4, -11.2), "AA/TI": (-0.7, -2.6),
    //     "TA/AI": (-1.3, -4.6), "CI/GA": (2.6, 8.9), "GI/CA": (-7.8, -21.1),
    //     "CA/GI": (-7.0, -20.0), "GA/CI": (-7.6, -20.2),
    //     "AI/TT": (0.49, -0.7), "TI/AT": (-6.5, -22.0), "AT/TI": (-5.6, -18.7),
    //     "TT/AI": (-0.8, -4.3), "CI/GT": (-1.0, -2.4), "GI/CT": (-3.5, -10.6),
    //     "CT/GI": (0.1, -1.0), "GT/CI": (-4.3, -12.1),
    //     "AI/TG": (-4.9, -15.8), "TI/AG": (-1.9, -8.5), "AG/TI": (0.1, -1.8),
    //     "TG/AI": (1.0, 1.0), "CI/GG": (7.1, 21.3), "GI/CG": (-1.1, -3.2),
    //     "CG/GI": (5.8, 16.9), "GG/CI": (-7.6, -22.0),
    //     "AI/TI": (-3.3, -11.9), "TI/AI": (0.1, -2.3), "CI/GI": (1.3, 3.0),
    //     "GI/CI": (-0.5, -1.3)}
}

/// Enthalpy (dH) and entropy (dS) tables based on dangling ends.
fn _dH_dS_de(nts: (Nucleotide, Nucleotide)) -> Option<(f32, f32)> {
    match nts {
        (A, A) => Some((0.2, 2.3)),
        (A, T) => Some((-7.2, -20.4)),
        (T, A) => Some((-7.2, -21.3)),
        (C, A) => Some((-8.5, -22.7)),
        (C, G) => Some((-10.6, -27.2)),
        (G, A) => Some((-8.2, -22.2)),
        (G, C) => Some((-9.8, -24.4)),
        (G, T) => Some((-8.4, -22.4)),
        (G, G) => Some((-8.0, -19.9)),
        _ => None,
    }

    // # Dangling ends table (DNA)
    // # Bommarito et al. (2000), Nucl Acids Res 28: 1929-1934
    // DNA_DE1 = {
    //     "AA/.T": (0.2, 2.3), "AC/.G": (-6.3, -17.1), "AG/.C": (-3.7, -10.0),
    //     "AT/.A": (-2.9, -7.6), "CA/.T": (0.6, 3.3), "CC/.G": (-4.4, -12.6),
    //     "CG/.C": (-4.0, -11.9), "CT/.A": (-4.1, -13.0), "GA/.T": (-1.1, -1.6),
    //     "GC/.G": (-5.1, -14.0), "GG/.C": (-3.9, -10.9), "GT/.A": (-4.2, -15.0),
    //     "TA/.T": (-6.9, -20.0), "TC/.G": (-4.0, -10.9), "TG/.C": (-4.9, -13.8),
    //     "TT/.A": (-0.2, -0.5),
    //     ".A/AT": (-0.7, -0.8), ".C/AG": (-2.1, -3.9), ".G/AC": (-5.9, -16.5),
    //     ".T/AA": (-0.5, -1.1), ".A/CT": (4.4, 14.9), ".C/CG": (-0.2, -0.1),
    //     ".G/CC": (-2.6, -7.4), ".T/CA": (4.7, 14.2), ".A/GT": (-1.6, -3.6),
    //     ".C/GG": (-3.9, -11.2), ".G/GC": (-3.2, -10.4), ".T/GA": (-4.1, -13.1),
    //     ".A/TT": (2.9, 10.4), ".C/TG": (-4.4, -13.1), ".G/TC": (-5.2, -15.0),
    //     ".T/TA": (-3.8, -12.6)}
}

/// Enthalpy (dH) and entropy (dS) based on nearest neighbors.
/// SantaLucia & Hicks, 2004, Table 1. Value are in kcal/Mol.
///
/// `neighbors` refers to the values between adjacent pairs of NTs.
/// To interpret this table, and come up with this result, search it in these
/// two manners:
/// A: Left to right, left of the slash
/// B: Right to left, right of the slash
/// You will find exactly one match using this approach.
fn dH_dS_neighbors(neighbors: (Nucleotide, Nucleotide)) -> (f32, f32) {
    match neighbors {
        (A, A) | (T, T) => (-7.6, -21.3),
        (A, T) => (-7.2, -20.4),
        (T, A) => (-7.2, -21.3),
        (C, A) | (T, G) => (-8.5, -22.7),
        (G, T) | (A, C) => (-8.4, -22.4),
        (C, T) | (A, G) => (-7.8, -21.0),
        (G, A) | (T, C) => (-8.2, -22.2),
        (C, G) => (-10.6, -27.2),
        (G, C) => (-9.8, -24.4),
        (G, G) | (C, C) => (-8.0, -19.9),
    }
}

/// Calculate a Tm correction term due to salt ions.
/// https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/MeltingTemp.py#L475
fn salt_correction(seq: &[Nucleotide], ion: &IonConcentrations) -> Option<f32> {
    // todo: Using Some casual defaults for now.
    // These are millimolar concentration of respective ions.
    let method = 5; // todo?

    let tris = 0.; // todo: Do we want this?

    if (5..=7).contains(&method) && seq.is_empty() {
        // return Err("sequence is missing (is needed to calculate GC content or sequence length).".into());
        return None;
    }

    let mut corr = 0.0;
    if method == 0 {
        return Some(corr);
    }

    // It appears that this section modifies the monovalent concentration with divalent values.
    let mut mon = ion.monovalent + tris / 2.0;
    let mon_molar = mon * 1e-3;
    let mg_molar = ion.divalent * 1e-3;

    if (ion.monovalent > 0.0 || ion.divalent > 0.0 || tris > 0.0 || ion.dntp > 0.0)
        && method != 7
        && ion.dntp < ion.divalent
    {
        mon += 120.0 * (ion.divalent - ion.dntp).sqrt();
    }

    if (1..=6).contains(&method) && mon_molar == 0.0 {
        // return Err("Total ion concentration of zero is not allowed in this method.".into());
        return None;
    }

    match method {
        1 => corr = 16.6 * mon_molar.log10(),
        2 => corr = 16.6 * (mon_molar / (1.0 + 0.7 * mon_molar)).log10(),
        3 => corr = 12.5 * mon_molar.log10(),
        4 => corr = 11.7 * mon_molar.log10(),
        5 => {
            corr = 0.368 * (seq.len() as f32 - 1.0) * mon_molar.ln();
        }
        6 => {
            let gc_fraction = calc_gc(seq);
            corr = ((4.29 * gc_fraction - 3.95) * 1e-5 * mon_molar.ln())
                + 9.40e-6 * mon_molar.ln().powi(2);
        }
        7 => {
            let (mut a, b, c, mut d, e, f, mut g) = (3.92, -0.911, 6.26, 1.42, -48.2, 52.5, 8.31);
            let mut mg_free = mg_molar;
            if ion.dntp > 0.0 {
                let dntps_molar = ion.dntp * 1e-3;
                let ka = 3e4;
                //  Free Mg2+ calculation
                mg_free = (-(ka * dntps_molar - ka * mg_free + 1.0)
                    + ((ka * dntps_molar - ka * mg_free + 1.0).powi(2) + 4.0 * ka * mg_free)
                        .sqrt())
                    / (2.0 * ka);
            }

            if mon > 0.0 {
                let r = (mg_free).sqrt() / mon_molar;
                if r < 0.22 {
                    let gc_fraction = calc_gc(seq);
                    corr = (4.29 * gc_fraction - 3.95) * 1e-5 * mon_molar.ln()
                        + 9.40e-6 * mon_molar.ln().powi(2);
                    return Some(corr);
                } else if r < 6.0 {
                    a = 3.92 * (0.843 - 0.352 * mon_molar.sqrt() * mon_molar.ln());
                    d = 1.42
                        * (1.279 - 4.03e-3 * mon_molar.ln() - 8.03e-3 * mon_molar.ln().powi(2));
                    g = 8.31 * (0.486 - 0.258 * mon_molar.ln() + 5.25e-3 * mon_molar.ln().powi(3));
                }

                let gc_fraction = calc_gc(seq);
                corr = (a
                    + b * mg_free.ln()
                    + gc_fraction * (c + d * mg_free.ln())
                    + (1.0 / (2.0 * (seq.len() as f32 - 1.0)))
                        * (e + f * mg_free.ln() + g * mg_free.ln().powi(2)))
                    * 1e-5;
            }
        }
        _ => return None,
    }

    Some(corr)
}

pub fn calc_tm(seq: &[Nucleotide], ion_concentrations: &IonConcentrations) -> Option<f32> {
    if seq.len() < MIN_PRIMER_LEN {
        return None;
    }

    // Inititial values. (S&H, Table 1)
    let mut dH = 0.2;
    let mut dS = -5.7;

    // If no GC content, apply additional values. (Table 1)
    if calc_gc(seq) < 0.001 {
        dH += 2.2;
        dS += 6.9;
    }

    // Add to dH and dS based on the terminal pair.
    {
        let term_pair = vec![seq[0], seq[seq.len() - 1]];
        let mut at_term_count = 0;
        // Constants for the CG term are 0, so we don't need it.
        for nt in term_pair {
            if nt == A || nt == T {
                at_term_count += 1;
            }
        }
        dH += 2.2 * at_term_count as f32;
        dS += 6.9 * at_term_count as f32;
    }

    for (i, nt) in seq.iter().enumerate() {
        if i + 1 >= seq.len() {
            break;
        }

        let neighbors = (*nt, seq[i + 1]);

        let (dH_nn, dS_nn) = dH_dS_neighbors(neighbors);
        dH += dH_nn;
        dS += dS_nn;
    }

    let C_T = ion_concentrations.primer * 1.0e-9;

    // println!("\n\ndH: {dH} dS: {dS}");

    if let Some(sc) = salt_correction(seq, ion_concentrations) {
        // Hard-coded for salt-correction method 5.
        dS += sc;
        // result += sc;

        // for saltcorr 6/7:
        // result = 1. / (1. / (result + 273.15) + sc) - 273.15
    } else {
        eprintln!("Error calculating salt correction.");
    }

    // SantaLucia and Hicks, Equation 3. Note the C_T / 2 vice / 4, due to double-stranded concentration.
    let result = (1_000. * dH) / (dS + R * (C_T / 2.).ln()) - 273.15;

    Some(result)
}
