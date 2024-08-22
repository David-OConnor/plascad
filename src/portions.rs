//! Code related to mixing portions. Used for performing quick mixing volume calculations.

use std::fmt::Display;

use bincode::{Decode, Encode};

#[derive(Default, Clone, Encode, Decode)]
pub struct PortionsState {
    pub solutions: Vec<Solution>,
}

#[derive(Default, Clone, Encode, Decode)]
pub struct Solution {
    pub name: String,
    /// Liters
    pub total_volume: f32,
    pub reagents: Vec<Reagent>,
    pub sub_solns: Vec<Solution>,
    /// Volatile; not to be added to directly.
    reagents_sub_solns: Vec<Reagent>,
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
                } else if *v > 0.001 {
                    format!("{:.2} mg", v * 1_000.)
                } else {
                    format!("{:.2} μg", v * 1_000_000.)
                }
            }
            // todo: Be careful about this calc being called frequently.
            Self::Volume(v) => {
                if *v > 1. {
                    format!("{:.2} L", v)
                } else if *v > 0.001 {
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

impl Reagent {
    pub fn calc_amount(&mut self, total_volume: f32) {
        // mol = mol/L x L
        let moles_req = self.molarity * total_volume;

        self.amount_calc = match self.prep {
            // g = g/mol x mol
            ReagentPrep::Mass => AmountCalculated::Mass(self.type_.weight() * moles_req),
            // L = mol / mol/L:
            ReagentPrep::Volume(reagent_molarity) => {
                AmountCalculated::Volume(moles_req / reagent_molarity)
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
    TrisHcl,
    Imidazole,
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
            Self::TrisHcl => 121.14,
            Self::Imidazole => 68.08,
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
            Self::TrisHcl => "Tris-HCl".to_owned(),
            Self::Imidazole => "Imidazole".to_owned(),
            Self::Custom(weight) => format!("Custom. Weight: {weight})"),
            Self::Solution(_) => format!("Solution").to_owned(), // todo
        };

        write!(f, "{}", str)
    }
}
