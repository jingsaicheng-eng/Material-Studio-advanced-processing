#!perl
use strict;
use warnings;
use MaterialsScript qw(:all);

# -------------------------------------------------------------------------
# Materials Studio script:
#   Calculate Born-charge polarization from a CASTEP .castep output.
#
# Usage in Materials Studio:
#   1. Put this script in the MS Script window.
#   2. Make sure BaTiO3_PhonDisp.castep is open/imported in the project, or
#      set CastepFilePath to the full file path.
#   3. Run the script.
#
# Output documents in the MS project:
#   <input>_polarization.txt
#   <input>_polarization.std
# -------------------------------------------------------------------------

# If the .castep file is already in the Project Explorer, use this name.
my $CastepDocumentName = "BaTiO3_PhonDisp.castep";

# Optional full path. Leave empty to use the Project Explorer document above.
my $CastepFilePath = "";

my $NUM_RE = qr/[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?/;

my ($input_name, @lines) = get_castep_lines($CastepDocumentName, $CastepFilePath);
my $base_name = base_from_name($input_name);

my ($a, $b, $c, $alpha, $beta, $gamma, $volume, $vectors_ref) =
    parse_lattice(\@lines);
my @atoms = parse_fractional_coordinates(\@lines);
my %serial_to_key;
for (my $i = 0; $i < @atoms; ++$i) {
    $serial_to_key{$i + 1} = $atoms[$i]->{key};
}
my %born = parse_born_charges(\@lines, \%serial_to_key);

die "Cannot find lattice parameters a/b/c in $input_name.\n"
    unless defined($a) && defined($b) && defined($c);
die "Cannot find or calculate cell volume in $input_name.\n"
    unless defined($volume) && $volume > 0.0;
die "Cannot find atom fractional coordinates in $input_name.\n"
    unless @atoms;
die "Cannot find Born effective charges in $input_name.\n"
    unless %born;

# P4mm-BaTiO3 centrosymmetric reference structure.
# Edit this table if your reference phase or atom numbering is different.
my %ref = (
    "O1"  => [ 0.000000, 0.000000, 0.500000 ],
    "O2"  => [ 0.000000, 0.500000, 0.000000 ],
    "O3"  => [ -0.500000, 0.000000, 0.000000 ],
    "Ti1" => [ 0.000000, 0.000000, 0.000000 ],
    "Ba1" => [ 0.500000, 0.500000, 0.500000 ],
);

my $text_doc = Documents->New(unique_document_name("${base_name}_polarization", "txt"));
my $std = Documents->New(unique_document_name("${base_name}_polarization", "std"));
setup_study_table($std);

my $factor = 1602.176634 / $volume;
my ($sum_x, $sum_y, $sum_z) = (0.0, 0.0, 0.0);

emit($text_doc, "Input CASTEP document: $input_name\n\n");
emit($text_doc, "Lattice parameters:\n");
emit($text_doc, sprintf("a = %.7f A\n", $a));
emit($text_doc, sprintf("b = %.7f A\n", $b));
emit($text_doc, sprintf("c = %.7f A\n", $c));
emit($text_doc, sprintf("alpha = %.7f deg\n", $alpha)) if defined $alpha;
emit($text_doc, sprintf("beta  = %.7f deg\n", $beta))  if defined $beta;
emit($text_doc, sprintf("gamma = %.7f deg\n", $gamma)) if defined $gamma;
emit($text_doc, sprintf("Volume = %.6f A^3\n", $volume));
emit($text_doc, sprintf("Unit conversion factor = 1602.176634 / Volume = %.8f\n\n", $factor));

emit($text_doc, "Atom-by-atom polarization contribution using full Born tensor:\n");
emit($text_doc, "Atom       dfx        dfy        dfz        dx_A       dy_A       dz_A       Cx_eA      Cy_eA      Cz_eA\n");

my $row = 0;
foreach my $atom (@atoms) {
    my $key = $atom->{key};

    die "No reference coordinate for '$key'. Edit %ref in the script.\n"
        unless exists $ref{$key};
    die "No Born effective charge for '$key'. Check CASTEP Born section.\n"
        unless exists $born{$key};

    my ($fx_ref, $fy_ref, $fz_ref) = @{$ref{$key}};
    my ($fx, $fy, $fz) = ($atom->{fx}, $atom->{fy}, $atom->{fz});

    my $dfx = min_image($fx - $fx_ref);
    my $dfy = min_image($fy - $fy_ref);
    my $dfz = min_image($fz - $fz_ref);

    my ($dx, $dy, $dz) = cartesian_displacement(
        $dfx, $dfy, $dfz, $a, $b, $c, $vectors_ref
    );

    my $z = $born{$key};
    my $cx = $z->[0][0] * $dx + $z->[0][1] * $dy + $z->[0][2] * $dz;
    my $cy = $z->[1][0] * $dx + $z->[1][1] * $dy + $z->[1][2] * $dz;
    my $cz = $z->[2][0] * $dx + $z->[2][1] * $dy + $z->[2][2] * $dz;

    $sum_x += $cx;
    $sum_y += $cy;
    $sum_z += $cz;

    emit($text_doc, sprintf(
        "%-4s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
        $key, $dfx, $dfy, $dfz, $dx, $dy, $dz, $cx, $cy, $cz
    ));

    write_atom_row($std, $row, $atom, $fx_ref, $fy_ref, $fz_ref,
                   $dfx, $dfy, $dfz, $dx, $dy, $dz, $z, $cx, $cy, $cz,
                   $factor);
    ++$row;
}

my $px = $factor * $sum_x;
my $py = $factor * $sum_y;
my $pz = $factor * $sum_z;

emit($text_doc, "\nSummation:\n");
emit($text_doc, sprintf("Sum Cx = %.8f e*A\n", $sum_x));
emit($text_doc, sprintf("Sum Cy = %.8f e*A\n", $sum_y));
emit($text_doc, sprintf("Sum Cz = %.8f e*A\n", $sum_z));

emit($text_doc, "\nPolarization:\n");
emit($text_doc, sprintf("Px = %.6f microC/cm^2\n", $px));
emit($text_doc, sprintf("Py = %.6f microC/cm^2\n", $py));
emit($text_doc, sprintf("Pz = %.6f microC/cm^2\n", $pz));
emit($text_doc, sprintf("P = (%.6f, %.6f, %.6f) microC/cm^2\n", $px, $py, $pz));

write_summary_rows($std, $row + 2, $a, $b, $c, $volume, $factor,
                   $sum_x, $sum_y, $sum_z, $px, $py, $pz);

$text_doc->Save;
$std->Save;

print "Finished Born-charge polarization analysis.\n";
print "Text report: " . $text_doc->Name . "\n";
print "Study table: " . $std->Name . "\n";

sub get_castep_lines {
    my ($doc_name, $file_path) = @_;
    my $doc = undef;

    if (defined($file_path) && $file_path ne "") {
        $doc = eval { Documents->Import($file_path); };
        die "Cannot import CASTEP file '$file_path'.\n$@\n" unless defined $doc;
    } else {
        $doc = eval { $Documents{$doc_name}; };
        if (!defined $doc) {
            $doc = eval { Documents->Import($doc_name); };
        }
    }

    die "Cannot find/import CASTEP document '$doc_name'. "
      . "Open it in Materials Studio or set CastepFilePath.\n"
        unless defined $doc;

    my @text_lines;
    my $ok = eval {
        foreach my $line (@{$doc->Lines}) {
            push @text_lines, $line;
        }
        1;
    };

    die "The document '" . $doc->Name . "' is not readable as a Text document.\n"
        unless $ok && @text_lines;

    return ($doc->Name, @text_lines);
}

sub parse_lattice {
    my ($lines) = @_;
    my ($a, $b, $c, $alpha, $beta, $gamma, $volume);
    my @vectors;

    for (my $i = 0; $i < @$lines; ++$i) {
        my $line = $lines->[$i];

        if ($line =~ /^\s*a\s*=\s*($NUM_RE)\b/i) {
            $a = to_num($1);
            $alpha = to_num($1) if $line =~ /\balpha\s*=?\s*($NUM_RE)/i;
        }
        if ($line =~ /^\s*b\s*=\s*($NUM_RE)\b/i) {
            $b = to_num($1);
            $beta = to_num($1) if $line =~ /\bbeta\s*=?\s*($NUM_RE)/i;
        }
        if ($line =~ /^\s*c\s*=\s*($NUM_RE)\b/i) {
            $c = to_num($1);
            $gamma = to_num($1) if $line =~ /\bgamma\s*=?\s*($NUM_RE)/i;
        }
        if ($line =~ /Current\s+cell\s+volume\s*=\s*($NUM_RE)/i) {
            $volume = to_num($1);
        }

        if ($line =~ /Real\s+Lattice/i) {
            my @tmp;
            for (my $j = $i + 1; $j < @$lines && @tmp < 3 && $j < $i + 10; ++$j) {
                my @n = numbers_from_line($lines->[$j]);
                next unless @n >= 3;
                push @tmp, [ @n[0, 1, 2] ];
            }
            if (@tmp == 3) {
                @vectors = @tmp;
                $a = vector_length($vectors[0]);
                $b = vector_length($vectors[1]);
                $c = vector_length($vectors[2]);
                $volume = abs(det3(\@vectors)) unless defined $volume;
            }
        }
    }

    if (!defined($volume) && defined($a) && defined($b) && defined($c)) {
        if (defined($alpha) && defined($beta) && defined($gamma)) {
            my $ca = cos(deg2rad($alpha));
            my $cb = cos(deg2rad($beta));
            my $cg = cos(deg2rad($gamma));
            $volume = $a * $b * $c
                    * sqrt(1.0 + 2.0 * $ca * $cb * $cg
                         - $ca * $ca - $cb * $cb - $cg * $cg);
        } else {
            $volume = $a * $b * $c;
        }
    }

    my $vectors_ref = @vectors == 3 ? \@vectors : undef;
    return ($a, $b, $c, $alpha, $beta, $gamma, $volume, $vectors_ref);
}

sub parse_fractional_coordinates {
    my ($lines) = @_;
    my @best_atoms;

    for (my $i = 0; $i < @$lines; ++$i) {
        next unless $lines->[$i] =~ /Fractional\s+coordinates\s+of\s+atoms/i;

        my @current;
        for (my $j = $i + 1; $j < @$lines; ++$j) {
            my $line = $lines->[$j];
            last if @current && $line =~ /Details\s+of\s+Species/i;

            if ($line =~ /^\s*x?\s*([A-Z][a-z]?)\s+(\d+)\s+($NUM_RE)\s+($NUM_RE)\s+($NUM_RE)(?:\s+x)?\s*$/) {
                push @current, atom_record($1, $2, $3, $4, $5);
                next;
            }
            if ($line =~ /^\s*([A-Z][a-z]?)\s+(\d+)\s+($NUM_RE)\s+($NUM_RE)\s+($NUM_RE)\s*$/) {
                push @current, atom_record($1, $2, $3, $4, $5);
                next;
            }
        }

        @best_atoms = @current if @current;
    }

    return @best_atoms;
}

sub parse_born_charges {
    my ($lines, $serial_to_key) = @_;
    my %born;
    my $in_born = 0;

    for (my $i = 0; $i < @$lines; ++$i) {
        my $line = $lines->[$i];
        if ($line =~ /^\s*Born\s+Effective\s+Charges\s*$/i) {
            $in_born = 1;
            next;
        }
        next unless $in_born;

        last if %born && $line =~ /^\s*=+\s*$/;
        last if %born && $line =~ /Symmetrised\s+Forces|Dielectric|Permittivity|Phonon|Raman|Infrared/i;

        if ($line =~ /^\s*x?\s*([A-Z][a-z]?)\s+(\d+)\s+($NUM_RE)\s+($NUM_RE)\s+($NUM_RE)/) {
            my $key = $1 . $2;
            my @rows = ([ to_num($3), to_num($4), to_num($5) ],
                        next_numeric_rows($lines, $i + 1, 2));
            $born{$key} = [ @rows[0, 1, 2] ] if @rows >= 3;
            next;
        }

        if ($line =~ /^\s*x?\s*([A-Z][a-z]?)\s+(\d+)\s*x?\s*$/) {
            my $key = $1 . $2;
            my @rows = next_numeric_rows($lines, $i + 1, 3);
            $born{$key} = [ @rows[0, 1, 2] ] if @rows >= 3;
            next;
        }

        if ($line =~ /^\s*(\d+)\s+([A-Z][a-z]?)\s+($NUM_RE)\s+($NUM_RE)\s+($NUM_RE)/) {
            my $serial = $1;
            my $key = exists $serial_to_key->{$serial}
                ? $serial_to_key->{$serial}
                : $2 . $serial;
            my @rows = ([ to_num($3), to_num($4), to_num($5) ],
                        next_numeric_rows($lines, $i + 1, 2));
            $born{$key} = [ @rows[0, 1, 2] ] if @rows >= 3;
            next;
        }
    }

    return %born;
}

sub atom_record {
    my ($elem, $num, $fx, $fy, $fz) = @_;
    return {
        key  => $elem . $num,
        elem => $elem,
        num  => $num,
        fx   => to_num($fx),
        fy   => to_num($fy),
        fz   => to_num($fz),
    };
}

sub next_numeric_rows {
    my ($lines, $start, $need) = @_;
    my @rows;

    for (my $j = $start; $j < @$lines && @rows < $need && $j < $start + 12; ++$j) {
        my @n = numbers_from_line($lines->[$j]);
        next unless @n >= 3;
        push @rows, [ @n[-3, -2, -1] ];
    }

    return @rows;
}

sub setup_study_table {
    my ($std) = @_;
    my @headings = qw(
        Atom Element Number
        fx_ref fy_ref fz_ref fx fy fz dfx dfy dfz dx_A dy_A dz_A
        Zxx Zxy Zxz Zyx Zyy Zyz Zzx Zzy Zzz
        Cx_eA Cy_eA Cz_eA Px_uC_cm2 Py_uC_cm2 Pz_uC_cm2
    );
    for (my $i = 0; $i < @headings; ++$i) {
        $std->ColumnHeading($i) = $headings[$i];
    }
}

sub write_atom_row {
    my ($std, $row, $atom, $fx_ref, $fy_ref, $fz_ref,
        $dfx, $dfy, $dfz, $dx, $dy, $dz, $z, $cx, $cy, $cz,
        $factor) = @_;

    my @values = (
        $atom->{key}, $atom->{elem}, $atom->{num},
        $fx_ref, $fy_ref, $fz_ref, $atom->{fx}, $atom->{fy}, $atom->{fz},
        $dfx, $dfy, $dfz, $dx, $dy, $dz,
        $z->[0][0], $z->[0][1], $z->[0][2],
        $z->[1][0], $z->[1][1], $z->[1][2],
        $z->[2][0], $z->[2][1], $z->[2][2],
        $cx, $cy, $cz,
        $factor * $cx, $factor * $cy, $factor * $cz
    );
    for (my $i = 0; $i < @values; ++$i) {
        $std->Cell($row, $i) = $values[$i];
    }
}

sub write_summary_rows {
    my ($std, $row, $a, $b, $c, $volume, $factor,
        $sum_x, $sum_y, $sum_z, $px, $py, $pz) = @_;

    $std->Cell($row, 0) = "Summary";
    $std->Cell($row + 1, 0) = "a_A";      $std->Cell($row + 1, 1) = $a;
    $std->Cell($row + 2, 0) = "b_A";      $std->Cell($row + 2, 1) = $b;
    $std->Cell($row + 3, 0) = "c_A";      $std->Cell($row + 3, 1) = $c;
    $std->Cell($row + 4, 0) = "Volume_A3";$std->Cell($row + 4, 1) = $volume;
    $std->Cell($row + 5, 0) = "Factor";   $std->Cell($row + 5, 1) = $factor;
    $std->Cell($row + 6, 0) = "Sum_Cx_eA";$std->Cell($row + 6, 1) = $sum_x;
    $std->Cell($row + 7, 0) = "Sum_Cy_eA";$std->Cell($row + 7, 1) = $sum_y;
    $std->Cell($row + 8, 0) = "Sum_Cz_eA";$std->Cell($row + 8, 1) = $sum_z;
    $std->Cell($row + 9, 0) = "Px_uC_cm2";$std->Cell($row + 9, 1) = $px;
    $std->Cell($row + 10, 0) = "Py_uC_cm2";$std->Cell($row + 10, 1) = $py;
    $std->Cell($row + 11, 0) = "Pz_uC_cm2";$std->Cell($row + 11, 1) = $pz;
}

sub numbers_from_line {
    my ($line) = @_;
    return map { to_num($_) } ($line =~ /($NUM_RE)/g);
}

sub to_num {
    my ($value) = @_;
    $value =~ tr/Dd/Ee/;
    return 0.0 + $value;
}

sub min_image {
    my ($df) = @_;
    $df -= 1.0 while $df > 0.5;
    $df += 1.0 while $df < -0.5;
    return $df;
}

sub cartesian_displacement {
    my ($dfx, $dfy, $dfz, $a, $b, $c, $vectors) = @_;
    if (defined $vectors && @$vectors == 3) {
        my $dx = $dfx * $vectors->[0][0] + $dfy * $vectors->[1][0] + $dfz * $vectors->[2][0];
        my $dy = $dfx * $vectors->[0][1] + $dfy * $vectors->[1][1] + $dfz * $vectors->[2][1];
        my $dz = $dfx * $vectors->[0][2] + $dfy * $vectors->[1][2] + $dfz * $vectors->[2][2];
        return ($dx, $dy, $dz);
    }
    return ($a * $dfx, $b * $dfy, $c * $dfz);
}

sub vector_length {
    my ($v) = @_;
    return sqrt($v->[0] * $v->[0] + $v->[1] * $v->[1] + $v->[2] * $v->[2]);
}

sub det3 {
    my ($m) = @_;
    return
        $m->[0][0] * ($m->[1][1] * $m->[2][2] - $m->[1][2] * $m->[2][1])
      - $m->[0][1] * ($m->[1][0] * $m->[2][2] - $m->[1][2] * $m->[2][0])
      + $m->[0][2] * ($m->[1][0] * $m->[2][1] - $m->[1][1] * $m->[2][0]);
}

sub deg2rad {
    my ($deg) = @_;
    return $deg * 3.14159265358979323846 / 180.0;
}

sub base_from_name {
    my ($name) = @_;
    $name =~ s/\\/\//g;
    $name =~ s!.*/!!;
    $name =~ s/\.[^.]+$//;
    $name =~ s/[^A-Za-z0-9_-]+/_/g;
    $name =~ s/^_+|_+$//g;
    return $name ne "" ? $name : "CASTEP";
}

sub unique_document_name {
    my ($base, $extension) = @_;
    $base =~ s/[^A-Za-z0-9_-]+/_/g;
    $base =~ s/^_+|_+$//g;
    $base = "MSDoc" if $base eq "";

    my $name = "$base.$extension";
    return $name unless document_exists($name);

    for (my $i = 1; $i < 10000; ++$i) {
        $name = sprintf("%s_%03d.%s", $base, $i, $extension);
        return $name unless document_exists($name);
    }
    die "Could not create a unique .$extension document name for $base.\n";
}

sub document_exists {
    my ($name) = @_;
    my $doc = eval { $Documents{$name}; };
    return defined($doc) ? 1 : 0;
}

sub emit {
    my ($text_doc, $text) = @_;
    print $text;
    $text_doc->Append($text);
}
