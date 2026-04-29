#!perl

use strict;
use warnings;
use Getopt::Long;
use MaterialsScript qw(:all);

# Planar average of a CASTEP field along X, Y, Z, AUTO, or ALL.
# X: average over YZ planes; Y: average over XZ planes; Z: average over XY planes.
# AUTO chooses the longest lattice direction, which is usually the vacuum/slab normal.
# Author Jingsai Cheng
# E-mail jingsaicheng@gmail.com
my $InputDocumentName = "1.xsd";
my $FieldName = "CASTEP SCF local potential";
my $Axis = "X";                  # X, Y, Z, AUTO, or ALL.
my $ProbeMode = "NearestGridPoint";
my $OutputStudyTableName = "";

# Coordinate ranges in Angstrom. Empty string means use the full cell range.
# For Axis=X, only Y/Z ranges restrict the average plane.
# For Axis=Y, only X/Z ranges restrict the average plane.
# For Axis=Z, only X/Y ranges restrict the average plane.
# The range along the profile axis trims which coordinates are written.
my $XMin_A = "";
my $XMax_A = "";
my $YMin_A = "4.95";
my $YMax_A = "8.1";
my $ZMin_A = "0";
my $ZMax_A = "25";

GetOptions(
    "doc=s" => \$InputDocumentName,
    "field=s" => \$FieldName,
    "axis=s" => \$Axis,
    "probe-mode=s" => \$ProbeMode,
    "out=s" => \$OutputStudyTableName,
    "xmin=f" => \$XMin_A,
    "xmax=f" => \$XMax_A,
    "ymin=f" => \$YMin_A,
    "ymax=f" => \$YMax_A,
    "zmin=f" => \$ZMin_A,
    "zmax=f" => \$ZMax_A,
);

my $doc = get_input_document($InputDocumentName);
my $field = get_named_field($doc, $FieldName);
my $field_probe = $field->CreateFieldProbe([ProbeMode => $ProbeMode]);

my @axis_labels = qw(X Y Z);
my @extents = (
    $field->GridExtent1,
    $field->GridExtent2,
    $field->GridExtent3,
);
my @lengths = lattice_lengths($doc);
my @range_min = ($XMin_A, $YMin_A, $ZMin_A);
my @range_max = ($XMax_A, $YMax_A, $ZMax_A);
normalize_ranges(\@range_min, \@range_max, \@lengths, \@axis_labels);
my @axes = axes_to_run($Axis, \@lengths);

my $std_name = $OutputStudyTableName;
$std_name = unique_document_name("FieldAvg", "std")
    if !defined($std_name) || $std_name eq "";
my $std = Documents->New($std_name);

$std->Cell(0, 0) = "Axis";
$std->Cell(0, 1) = "Index";
$std->Cell(0, 2) = "Coordinate_A";
$std->Cell(0, 3) = "AverageField";
$std->Cell(0, 4) = "SampleCount";
$std->Cell(0, 5) = "XRange_A";
$std->Cell(0, 6) = "YRange_A";
$std->Cell(0, 7) = "ZRange_A";
$std->Cell(0, 8) = "FieldName";
$std->Cell(0, 9) = "ProbeMode";

my $row = 1;
my $ok = eval {
    $field_probe->HasFieldLock = "Yes";

    foreach my $axis (@axes) {
        my $axis_index = axis_index($axis);
        print "\n# Axis $axis: coordinate_A average_field\n";
        $row = write_axis_average(
            $std, $row, $field_probe, $FieldName, $ProbeMode,
            $axis_index, \@extents, \@lengths, \@axis_labels,
            \@range_min, \@range_max
        );
    }

    1;
};
my $error = $@;
eval { $field_probe->HasFieldLock = "No"; };
die $error unless $ok;

$std->Save;
print "\nSaved planar average table: " . $std->Name . "\n";

sub get_input_document {
    my ($name) = @_;

    if (defined($name) && $name ne "") {
        my $document = eval { $Documents{$name}; };
        die "Could not find document '$name'. Check the XSD name in the project.\n"
            unless defined $document;
        return $document;
    }

    my $active = eval { Documents->ActiveDocument; };
    die "No active document. Set InputDocumentName or use -doc.\n"
        unless defined $active;
    return $active;
}

sub get_named_field {
    my ($document, $field_name) = @_;

    my $field_obj = eval { $document->AsymmetricUnit->Fields($field_name); };
    $field_obj = eval { $document->UnitCell->Fields($field_name); }
        unless defined $field_obj;
    $field_obj = eval { $document->Fields($field_name); }
        unless defined $field_obj;

    die "Could not find field '$field_name' in " . $document->Name . ".\n"
        unless defined $field_obj;
    return $field_obj;
}

sub lattice_lengths {
    my ($document) = @_;

    my $symmetry = eval { $document->SymmetryDefinition; };
    die "The document has no SymmetryDefinition; lattice lengths are required.\n"
        unless defined $symmetry;

    my @values = (
        eval { $symmetry->LengthA; },
        eval { $symmetry->LengthB; },
        eval { $symmetry->LengthC; },
    );

    foreach my $value (@values) {
        die "Could not read lattice LengthA/LengthB/LengthC.\n"
            unless defined($value) && $value > 0.0;
    }

    return @values;
}

sub normalize_ranges {
    my ($range_min, $range_max, $lengths, $axis_labels) = @_;

    for (my $i = 0; $i < 3; ++$i) {
        if (has_bound($range_min->[$i]) && has_bound($range_max->[$i])
            && $range_min->[$i] > $range_max->[$i]) {
            my $tmp = $range_min->[$i];
            $range_min->[$i] = $range_max->[$i];
            $range_max->[$i] = $tmp;
        }

        if (has_bound($range_min->[$i]) && $range_min->[$i] < 0.0) {
            die "$axis_labels->[$i]Min_A must be >= 0.0 A.\n";
        }

        if (has_bound($range_max->[$i]) && $range_max->[$i] > $lengths->[$i]) {
            die "$axis_labels->[$i]Max_A must be <= lattice length "
                . $lengths->[$i] . " A.\n";
        }
    }
}

sub has_bound {
    my ($value) = @_;
    return defined($value) && $value ne "";
}

sub axes_to_run {
    my ($axis_text, $lengths) = @_;
    $axis_text = uc($axis_text || "Z");

    return qw(X Y Z) if $axis_text eq "ALL";

    if ($axis_text eq "AUTO") {
        my $best = 0;
        for (my $i = 1; $i < 3; ++$i) {
            $best = $i if $lengths->[$i] > $lengths->[$best];
        }
        print "AUTO axis selected: " . (qw(X Y Z))[$best] . "\n";
        return ((qw(X Y Z))[$best]);
    }

    die "Axis must be X, Y, Z, AUTO, or ALL.\n"
        unless $axis_text =~ /^[XYZ]$/;
    return ($axis_text);
}

sub axis_index {
    my ($axis) = @_;
    return 0 if $axis eq "X";
    return 1 if $axis eq "Y";
    return 2 if $axis eq "Z";
    die "Unknown axis '$axis'.\n";
}

sub write_axis_average {
    my ($std, $row, $probe, $field_name, $probe_mode,
        $axis_index, $extents, $lengths, $axis_labels,
        $range_min, $range_max) = @_;

    my @plane_indices = grep { $_ != $axis_index } (0, 1, 2);
    my $axis_extent = $extents->[$axis_index];
    my $plane_a_extent = $extents->[$plane_indices[0]];
    my $plane_b_extent = $extents->[$plane_indices[1]];
    my $axis_label = $axis_labels->[$axis_index];

    die "Field grid extent is zero along $axis_label.\n"
        unless $axis_extent > 0 && $plane_a_extent > 0 && $plane_b_extent > 0;

    for (my $i = 0; $i < $axis_extent; ++$i) {
        my $coordinate =
            voxel_coordinate($i, $axis_extent, $lengths->[$axis_index]);
        next unless coordinate_in_range(
            $coordinate,
            $range_min->[$axis_index],
            $range_max->[$axis_index]
        );

        my $total = 0.0;
        my $count = 0;

        for (my $j = 0; $j < $plane_a_extent; ++$j) {
            my $plane_a_index = $plane_indices[0];
            my $plane_a_coord =
                voxel_coordinate($j, $extents->[$plane_a_index],
                                 $lengths->[$plane_a_index]);
            next unless coordinate_in_range(
                $plane_a_coord,
                $range_min->[$plane_a_index],
                $range_max->[$plane_a_index]
            );

            for (my $k = 0; $k < $plane_b_extent; ++$k) {
                my $plane_b_index = $plane_indices[1];
                my $plane_b_coord =
                    voxel_coordinate($k, $extents->[$plane_b_index],
                                     $lengths->[$plane_b_index]);
                next unless coordinate_in_range(
                    $plane_b_coord,
                    $range_min->[$plane_b_index],
                    $range_max->[$plane_b_index]
                );

                my @voxel = (0, 0, 0);
                $voxel[$axis_index] = $i;
                $voxel[$plane_indices[0]] = $j;
                $voxel[$plane_indices[1]] = $k;

                $probe->ProbeVoxelPosition = Point(
                    X => $voxel[0],
                    Y => $voxel[1],
                    Z => $voxel[2],
                );

                my $value = eval { $probe->FieldValue; };
                next if $@ || !defined $value;
                $total += $value;
                ++$count;
            }
        }

        my $average = $count > 0 ? $total / $count : "";

        print "$coordinate  $average\n";
        $std->Cell($row, 0) = $axis_label;
        $std->Cell($row, 1) = $i;
        $std->Cell($row, 2) = $coordinate;
        $std->Cell($row, 3) = $average;
        $std->Cell($row, 4) = $count;
        $std->Cell($row, 5) = range_text($range_min, $range_max, $lengths, 0);
        $std->Cell($row, 6) = range_text($range_min, $range_max, $lengths, 1);
        $std->Cell($row, 7) = range_text($range_min, $range_max, $lengths, 2);
        $std->Cell($row, 8) = $field_name;
        $std->Cell($row, 9) = $probe_mode;
        ++$row;
    }

    return $row;
}

sub voxel_coordinate {
    my ($index, $extent, $length) = @_;
    return 0.0 unless defined($extent) && $extent > 0;
    return $index * $length / $extent;
}

sub coordinate_in_range {
    my ($coordinate, $min, $max) = @_;
    return 0 if has_bound($min) && $coordinate < $min;
    return 0 if has_bound($max) && $coordinate > $max;
    return 1;
}

sub range_text {
    my ($range_min, $range_max, $lengths, $axis_index) = @_;
    my $min = has_bound($range_min->[$axis_index])
        ? $range_min->[$axis_index]
        : 0.0;
    my $max = has_bound($range_max->[$axis_index])
        ? $range_max->[$axis_index]
        : $lengths->[$axis_index];
    return "$min:$max";
}

sub unique_document_name {
    my ($prefix, $extension) = @_;
    $prefix =~ s/[^A-Za-z0-9]+//g;
    $prefix = "D" if $prefix eq "";

    for (my $attempt = 0; $attempt < 1000; ++$attempt) {
        my $suffix = $attempt == 0 ? "" : $attempt;
        my $base = substr($prefix . $suffix, 0, 24);
        my $name = "$base.$extension";
        my $existing = eval { $Documents{$name}; };
        return $name unless defined $existing;
    }

    die "Could not create a unique .$extension document name.\n";
}
