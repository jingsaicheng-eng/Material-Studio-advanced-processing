#!perl

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

# Author Jingsai Cheng
# E-mail jingsaicheng@gmail.com
# Change CASTEP cdd fiel to VASP CHGDIFF
# --- 1. Get current document and grid field data ---
my $doc = $Documents{"BN-P.xsd"};
my $field = $doc->AsymmetricUnit->Fields("CASTEP density difference from sets");
my $fieldProbe = $field->CreateFieldProbe([ProbeMode => "NearestGridPoint"]);

my $lattice = $doc->SymmetryDefinition;
my $xExtent = $field->GridExtent1;
my $yExtent = $field->GridExtent2;
my $zExtent = $field->GridExtent3;
my $numAtoms = $doc->UnitCell->Atoms->Count;

# --- 2. Create document object using MS native API ---
my $text = Documents->New("BN-P_diff.txt"); 

# --- 3. Write Cube file header ---
$text->Append(sprintf "Generated from Materials Studio Perl Script\n");
$text->Append(sprintf "Charge Density Difference Field\n");
$text->Append(sprintf "%5d %12.6f %12.6f %12.6f\n", $numAtoms, 0, 0, 0);

my $vecA = $lattice->VectorA;
my $vecB = $lattice->VectorB;
my $vecC = $lattice->VectorC;

# Write lattice vectors (Grid dimensions and step size, unit: Bohr)
$text->Append(sprintf "%5d %12.6f %12.6f %12.6f\n", $xExtent, $vecA->X*1.88973/$xExtent, $vecA->Y*1.88973/$xExtent, $vecA->Z*1.88973/$xExtent);
$text->Append(sprintf "%5d %12.6f %12.6f %12.6f\n", $yExtent, $vecB->X*1.88973/$yExtent, $vecB->Y*1.88973/$yExtent, $vecB->Z*1.88973/$yExtent);
$text->Append(sprintf "%5d %12.6f %12.6f %12.6f\n", $zExtent, $vecC->X*1.88973/$zExtent, $vecC->Y*1.88973/$zExtent, $vecC->Z*1.88973/$zExtent);

# --- 4. Write atom coordinates ---
foreach my $atom (@{$doc->UnitCell->Atoms}) {
    $text->Append(sprintf "%5d %12.6f %12.6f %12.6f %12.6f\n", 
        $atom->AtomicNumber, 0, $atom->X*1.88973, $atom->Y*1.88973, $atom->Z*1.88973);
}

# --- 5. Loop over 3D grid and write density data ---
$fieldProbe->HasFieldLock = "Yes";
my $count = 0;
my $line_buffer = "";

for (my $i=0; $i<$xExtent; ++$i) {
    for (my $j=0; $j<$yExtent; ++$j) {
        for (my $k=0; $k<$zExtent; ++$k) {
            $fieldProbe->ProbeVoxelPosition = Point(X => $i, Y => $j, Z => $k);
            
            # Concatenate data for the current line
            $line_buffer .= sprintf(" %13.5E", $fieldProbe->FieldValue);
            $count++;
            
            # Cube format requirement: Max 6 data points per line
            if ($count % 6 == 0) { 
                $text->Append(sprintf "%s\n", $line_buffer);
                $line_buffer = "";
            }
        }
        # Force newline if current layer ends with less than 6 data points
        if ($count % 6 != 0) { 
            $text->Append(sprintf "%s\n", $line_buffer);
            $line_buffer = ""; 
            $count = 0; 
        }
    }
}
$fieldProbe->HasFieldLock = "No";

# Note: Remember to press Ctrl+S to save the generated file after execution.
