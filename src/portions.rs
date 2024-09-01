//! Code related to mixing portions. Used for performing quick mixing volume calculations.

// todo: Consider including pKa info, and including tools to balance pH.

use std::fmt::Display;

use bincode::{Decode, Encode};

#[derive(Default, Clone, Encode, Decode)]
pub struct PortionsState {
    pub solutions: Vec<Solution>,
    pub media_input: MediaPrepInput,
    pub media_result: MediaPrep,
}

#[derive(Default, Clone, Encode, Decode)]
pub struct Solution {
    pub name: String,
    /// Liters
    pub total_volume: f32,
    pub reagents: Vec<Reagent>,
    pub sub_solns: Vec<Solution>,
    /// Volatile; not to be added to directly.
    pub reagents_sub_solns: Vec<Reagent>,
}

impl Solution {
    /// Find required amounts (mass or volume) for each reagent.
    pub fn calc_amounts(&mut self) {
        self.reagents_sub_solns = Vec::new();
        for sub_sol in &self.sub_solns {
            for reagent in &sub_sol.reagents {
                self.reagents_sub_solns.push(reagent.clone());
            }
        }

        for reagent in &mut self.reagents {
            reagent.calc_amount(self.total_volume);
        }
    }
}

#[derive(Clone, Copy, Encode, Decode)]
pub enum AmountCalculated {
    /// grams
    Mass(f32),
    /// Liters
    Volume(f32),
}

impl Display for AmountCalculated {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            Self::Mass(v) => {
                if *v > 1. {
                    format!("{:.2} g", v)
                } else if *v >= 0.001 {
                    format!("{:.2} mg", v * 1_000.)
                } else {
                    format!("{:.2} μg", v * 1_000_000.)
                }
            }
            // todo: Be careful about this calc being called frequently.
            Self::Volume(v) => {
                if *v > 1. {
                    format!("{:.2} L", v)
                } else if *v >= 0.001 {
                    format!("{:.2} mL", v * 1_000.)
                } else {
                    // todo: If you have precision issues, make your base unit mL.
                    format!("{:.2} μL", v * 1_000_000.)
                }
            }
        };

        write!(f, "{}", str)
    }
}

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum ReagentPrep {
    Mass,
    /// Inner: Molarity
    Volume(f32),
}

impl Display for ReagentPrep {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            Self::Mass => "Mass".to_owned(),
            // Self::Volume(molarity) => format!("Volume. Molarity: {molarity})"),
            Self::Volume(_molarity) => "Volume".to_owned(),
        };
        write!(f, "{}", str)
    }
}

#[derive(Clone, Encode, Decode)]
pub struct Reagent {
    pub type_: ReagentType,
    pub prep: ReagentPrep,
    /// Target moles per liter in the solution
    pub molarity: f32,
    /// Calculated result
    pub amount_calc: AmountCalculated,
}

impl Default for Reagent {
    fn default() -> Self {
        Self {
            type_: ReagentType::Custom(0.),
            prep: ReagentPrep::Mass,
            molarity: 0.,
            amount_calc: AmountCalculated::Mass(0.),
        }
    }
}

impl Reagent {
    pub fn calc_amount(&mut self, total_volume: f32) {
        // mol = mol/L x L
        let moles_req = self.molarity * total_volume;

        self.amount_calc = match self.prep {
            // g = g/mol x mol
            ReagentPrep::Mass => AmountCalculated::Mass(self.type_.weight() * moles_req),
            // L = mol / mol/L:
            ReagentPrep::Volume(reagent_molarity) => {
                if reagent_molarity.abs() < 0.00000001 {
                    AmountCalculated::Volume(0.)
                } else {
                    AmountCalculated::Volume(moles_req / reagent_molarity)
                }
            }
        };
    }
}

/// A collection of common reagents, where we have molecular weights.
#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum ReagentType {
    SodiumChloride,
    SodiumPhosphateMonobasic,
    SodiumPhosphateDibasic,
    SodiumPhosphateDibasicHeptahydrate,
    PotassiumPhosphateMonobasic,
    PotassiumPhosphateDibasic,
    TrisHcl,
    Iptg,
    Imidazole,
    Lysozyme,
    Mes,
    Bes,
    Tes,
    CitricAcid,
    Edta,
    HydrochloricAcid,
    // AceticAcid,
    SodiumHydroxide,
    Custom(f32), // Inner: Molecular weight
    /// Index of the solution in state.
    Solution(usize),
}

impl ReagentType {
    /// g/mol
    pub fn weight(&self) -> f32 {
        match self {
            Self::SodiumChloride => 58.44,
            Self::SodiumPhosphateMonobasic => 119.98,
            Self::SodiumPhosphateDibasic => 141.96,
            Self::SodiumPhosphateDibasicHeptahydrate => 268.10,
            Self::PotassiumPhosphateMonobasic => 136.08,
            Self::PotassiumPhosphateDibasic => 174.17,
            Self::TrisHcl => 121.14,
            Self::Iptg => 238.298,
            Self::Imidazole => 68.08,
            Self::Lysozyme => 14_388.,
            Self::Mes => 195.24,
            Self::Bes => 213.25,
            Self::Tes => 229.25,
            Self::CitricAcid => 192.12,
            Self::Edta => 292.24,
            Self::HydrochloricAcid => 36.46,
            Self::SodiumHydroxide => 40.,
            Self::Custom(weight) => *weight,
            Self::Solution(_) => 0., // todo?
        }
    }
}

impl Display for ReagentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            Self::SodiumChloride => "NaCl".to_owned(),
            Self::SodiumPhosphateMonobasic => "NaH₂PO₄ (Mono)".to_owned(),
            Self::SodiumPhosphateDibasic => "Na₂HPO₄ (Di)".to_owned(),
            Self::SodiumPhosphateDibasicHeptahydrate => "Na₂HPO₄·7H₂O".to_owned(),
            Self::PotassiumPhosphateMonobasic => "H₂KO₄P".to_owned(),
            Self::PotassiumPhosphateDibasic => "HK₂O₄P".to_owned(),
            Self::TrisHcl => "Tris-HCl".to_owned(),
            Self::Iptg => "IPTG".to_owned(),
            Self::Imidazole => "Imidazole".to_owned(),
            Self::Lysozyme => "Lysozyme".to_owned(),
            Self::Mes => "MES".to_owned(),
            Self::Bes => "BES".to_owned(),
            Self::Tes => "BES".to_owned(),
            Self::CitricAcid => "Citric acid".to_owned(),
            Self::Edta => "EDTA".to_owned(),
            Self::HydrochloricAcid => "HCl".to_owned(),
            Self::SodiumHydroxide => "NaOH".to_owned(),
            Self::Custom(_) => "Custom".to_owned(),
            Self::Solution(_) => format!("Solution").to_owned(), // todo
        };

        write!(f, "{}", str)
    }
}

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum PlateSize {
    /// 60mm diameter
    D60,
    D90,
    D100,
    D150,
}

impl PlateSize {
    /// Nominal amount of liquid. Note that these go with the square of the baseline. We use 7mL for 60mm
    /// plates as the baseline.
    pub fn volume(&self) -> f32 {
        match self {
            Self::D60 => 0.007,
            Self::D90 => 0.01575,
            Self::D100 => 0.01944,
            Self::D150 => 0.04375,
        }
    }
}

impl Display for PlateSize {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            Self::D60 => "60mm",
            Self::D90 => "90mm",
            Self::D100 => "100mm",
            Self::D150 => "150mm",
        };

        write!(f, "{}", str)
    }
}

#[derive(Clone, PartialEq, Encode, Decode)]
pub enum MediaPrepInput {
    Plates((PlateSize, usize)), // number of plates,
    Liquid(f32),                // volume
}

impl Display for MediaPrepInput {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            Self::Plates(_) => "Plates",
            Self::Liquid(_) => "Liquid culture",
        };

        write!(f, "{}", str)
    }
}

impl Default for MediaPrepInput {
    fn default() -> Self {
        Self::Plates((PlateSize::D90, 6))
    }
}

#[derive(Clone, Default, Encode, Decode)]
pub struct MediaPrep {
    pub water: f32,      // L
    pub food: f32,       // g
    pub agar: f32,       // g
    pub antibiotic: f32, // mL
}

/// Returns volume of water, grams of LB, grams of agar, mL of 1000x antibiotics.
pub fn media_prep(input: &MediaPrepInput) -> MediaPrep {
    let (volume, agar) = match input {
        MediaPrepInput::Plates((plate_size, num)) => {
            let volume = plate_size.volume() * *num as f32;
            (volume, volume * 15.)
        }
        MediaPrepInput::Liquid(volume) => (*volume, 0.),
    };

    MediaPrep {
        water: volume,
        food: volume * 25.,
        agar,
        antibiotic: volume,
    }
}
