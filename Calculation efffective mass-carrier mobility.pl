#!perl
use strict;
use warnings;
use Cwd qw(getcwd);
use File::Find;
use MaterialsScript qw(:all);

# -------------------------------------------------------------------------
# Linked calculation stage only.
#
# This script submits/generates all CASTEP raw data needed for:
#   effective-mass tensor + A/B/C deformation potentials + elastic moduli
#   + carrier mobility.
#
# It does not fit or calculate mobility.  Run the analysis stage afterward:
#   castep_mobility_linked_analysis.pl
# -------------------------------------------------------------------------

# ================= User input =================

my $InputDocumentName = "Si.xsd";
my $TemperatureK = 300.0;

my @Directions = ("A", "B", "C");
my @Strains = (-0.010, -0.005, 0.000, 0.005, 0.010);

my $RunInitialGeometryOptimization = "Yes";
my $InitialOptSeedSuffix = "OPT";
my $BandsNameKeyword = "BandStr";

my %CASTEPCommonSettings = (
    Quality => "Fine",
    UseCustomEnergyCutoff => "Yes",
    EnergyCutoff => 550,
    SCFConvergence => 2.0e-5,
    MaximumSCFCycles => 200,
    #UseDFTD => "Yes",
    KPointDerivation => "CustomGrid",
    ParameterA => 5,
    ParameterB => 5,
    ParameterC => 5,
    UseInsulatorDerivation => "Yes",
    Pseudopotentials => "Ultrasoft",
    XCFunctional => "PBE",
    #SpinTreatment => "Collinear",
    #UseDFTU => "Yes"
);

my %CASTEPGeometryOptimizationSettings = (
    %CASTEPCommonSettings,
    CellOptimization => "Full",
    OptimizationAlgorithm => "LBFGS"
);

my %CASTEPBandStructureSettings = (
    %CASTEPCommonSettings,
    CalculateBandStructure => "Dispersion",
    PropertiesKPointQuality => "Fine",
    BandStructureEmax => 10,
    BandStructureNumExtraBands => 12,
    BandStructureEnergyTolerance => 1.0e-5,
    BandStructureLineStyle => "Line"
);

my @BandsSearchDirs = (
    ".",
    getcwd(),
    "D:/msproject/bb_Files/Documents",
    "C:/Users/22903/Documents/New project",
    "C:/Users/22903/Desktop"
);
my @MassTensorSearchDirs = @BandsSearchDirs;

my $ReferenceSeedPrefix = "SiREF";
my $MassTensorSeedPrefix = "SiM3D";
my $MassGridHalfWidth = 2;
my $MassGridStepFrac = 0.010;
my $MassSCFKPointMPSpacing = "0.04 1/Ang";
my $MassCutOffEnergy = $CASTEPCommonSettings{EnergyCutoff} . " eV";
my $MassXCFunctional = $CASTEPCommonSettings{XCFunctional};
my $RunMassTensorFilesTask = "Yes";

my $MakeP1BeforeStrain = "Yes";
my $MakeP1BeforeFiles = "Yes";
my $SaveAllAtEnd = "Yes";
my $OutputStdName = "Linked_Mobility_CalcData.std";

my $FermiTolerance_eV = 1.0e-6;

# ================= Constants =================

my $HARTREE_TO_EV = 27.211386245988;
my $KCALMOL_TO_EV_PER_CELL = 1.0 / 23.060548867;

# ================= Main =================

Documents->SaveAllAtEnd = $SaveAllAtEnd;

my $source_doc = $Documents{$InputDocumentName};
die "Could not find Documents{\"$InputDocumentName\"}. Check the document name.\n"
    unless defined $source_doc;

$source_doc = optimize_initial_structure($source_doc);

my $lat0 = $source_doc->Lattice3D;
die "$InputDocumentName is not a 3D periodic document.\n" unless defined $lat0;

my $base_seed = strip_extension($InputDocumentName);
my ($a0, $b0, $c0, $v0) = ($lat0->LengthA, $lat0->LengthB, $lat0->LengthC, $lat0->CellVolume);

print "Linked calculation stage\n";
print "Input document: $InputDocumentName\n";
print "Working document: " . $source_doc->Name . ".xsd\n";
printf "Working cell: a=%g A, b=%g A, c=%g A, V=%g A^3\n\n", $a0, $b0, $c0, $v0;

my %mass_info = run_mass_tensor_calculation();
my @deformation_records = run_deformation_calculations();

write_linked_calc_std(\%mass_info, \@deformation_records);

print "\nLinked calculation data saved: $OutputStdName\n";
print "Next run: castep_mobility_linked_analysis.pl\n";

# ================= CASTEP calculations =================

sub optimize_initial_structure {
    my ($doc) = @_;
    return $doc unless $RunInitialGeometryOptimization =~ /^y/i;

    my $seed = optimized_seed_from_doc($doc);
    print "Running initial CASTEP geometry optimization: $seed.xsd, CellOptimization=$CASTEPGeometryOptimizationSettings{CellOptimization}\n";

    my $opt_doc = $doc->SaveAs("$seed.xsd");
    my $results = Modules->CASTEP->GeometryOptimization->Run($opt_doc, Settings(%CASTEPGeometryOptimizationSettings));

    my $optimized_doc = $opt_doc;
    eval {
        my $structure = $results->Structure;
        $optimized_doc = $structure if defined $structure;
    };
    $optimized_doc->Save;
    print "Initial geometry optimization finished. Using " . $optimized_doc->Name . ".xsd\n\n";
    return $optimized_doc;
}

sub run_castep_band_structure {
    my ($doc) = @_;
    return Modules->CASTEP->Energy->Run($doc, Settings(%CASTEPBandStructureSettings));
}

sub run_mass_tensor_calculation {
    my $ref_seed = short_seed($ReferenceSeedPrefix);
    my $ref_doc = $source_doc->SaveAs("$ref_seed.xsd");
    my $actual_ref_seed = $ref_doc->Name;

    print "Running reference BandStru calculation for $actual_ref_seed...\n";
    run_castep_band_structure($ref_doc);

    my $ref_bands = find_band_file($actual_ref_seed, $BandsNameKeyword, \@BandsSearchDirs);
    $ref_bands = find_band_file($ref_seed, $BandsNameKeyword, \@BandsSearchDirs) unless defined $ref_bands;
    die "Reference BandStru .bands was not found after CASTEP run.\n" unless defined $ref_bands;

    my $ref_data = read_castep_bands($ref_bands);
    my ($vbm, $cbm) = find_band_edges($ref_data);
    die "Could not locate VBM/CBM from reference BandStru .bands.\n" unless defined $vbm && defined $cbm;

    printf "Reference VBM: %.8f eV, kpoint %d, band %d\n",
        $vbm->{energy_eV}, $vbm->{kindex} + 1, $vbm->{band} + 1;
    printf "Reference CBM: %.8f eV, kpoint %d, band %d\n",
        $cbm->{energy_eV}, $cbm->{kindex} + 1, $cbm->{band} + 1;

    my @kpoints = build_mass_tensor_kpoint_grid($ref_data, $vbm, $cbm);
    my $kblock = spectral_kpoint_list_block(\@kpoints);

    my $mass_seed = short_seed($MassTensorSeedPrefix);
    my $mass_doc = $source_doc->SaveAs("$mass_seed.xsd");
    $mass_doc->MakeP1 if $MakeP1BeforeFiles =~ /^y/i;
    $mass_doc->Save;
    my $actual_mass_seed = $mass_doc->Name;

    my $cell_text = make_mass_tensor_cell_text($mass_doc, $kblock);
    my $param_text = make_mass_tensor_param_text();

    write_text_document("$actual_mass_seed.cell", $cell_text);
    write_text_document("$actual_mass_seed.param", $param_text);
    write_text_document("${actual_mass_seed}_SPECTRAL_KPOINT_LIST.txt", $kblock);

    my $project_dir = find_project_document_folder($actual_mass_seed, \@BandsSearchDirs);
    $project_dir = dirname($ref_bands) unless defined $project_dir;
    if (defined $project_dir) {
        write_external_text_file("$project_dir/$actual_mass_seed.cell", $cell_text);
        write_external_text_file("$project_dir/$actual_mass_seed.param", $param_text);
        write_external_text_file("$project_dir/${actual_mass_seed}_SPECTRAL_KPOINT_LIST.txt", $kblock);
        print "Wrote CASTEP Files input directly to project folder: $project_dir\n";
    }

    my $mass_bands = "";
    if ($RunMassTensorFilesTask =~ /^y/i) {
        print "Running CASTEP Files task for effective-mass tensor: $actual_mass_seed...\n";
        my $files_run_started = time;
        eval { Modules->CASTEP->Files->Run($mass_doc); };
        if ($@) {
            print "CASTEP Files task failed. The .cell/.param files were generated.\n";
            print "Error was: $@\n";
        } else {
            my $found = find_band_file($actual_mass_seed, "", \@MassTensorSearchDirs, $files_run_started - 2);
            $mass_bands = $found if defined $found;
        }
    }

    print "Mass tensor .bands: " . ($mass_bands ne "" ? $mass_bands : "<not found; fill it manually in $OutputStdName>") . "\n\n";

    return (
        reference_bands => $ref_bands,
        mass_bands => $mass_bands,
        mass_seed => $actual_mass_seed,
        generated_kpoints => scalar(@kpoints),
        vbm_eV => $vbm->{energy_eV},
        vbm_kpoint => $vbm->{kindex} + 1,
        vbm_band => $vbm->{band} + 1,
        cbm_eV => $cbm->{energy_eV},
        cbm_kpoint => $cbm->{kindex} + 1,
        cbm_band => $cbm->{band} + 1
    );
}

sub run_deformation_calculations {
    my @records;

    foreach my $dir (@Directions) {
        foreach my $strain (@Strains) {
            my $tag = strain_tag($strain);
            my $seed = short_deformation_seed($base_seed, $dir, $tag);

            print "--- $dir strain $strain ($tag) ---\n";
            my $work_doc = make_strained_copy($source_doc, "$seed.xsd", $dir, $strain);
            my $actual_seed = $work_doc->Name;
            print "Running CASTEP BandStru for $actual_seed...\n";

            my $results = run_castep_band_structure($work_doc);
            my ($total_energy_kcalmol, $total_energy_eV) = result_total_energy($results);

            my $bands_file = find_band_file($actual_seed, $BandsNameKeyword, \@BandsSearchDirs);
            $bands_file = find_band_file($seed, $BandsNameKeyword, \@BandsSearchDirs) unless defined $bands_file;
            die "No BandStru .bands found for $actual_seed. Add the project Documents folder to \@BandsSearchDirs.\n"
                unless defined $bands_file;

            my ($a, $b, $c, $v) = expected_cell_after_strain($dir, $strain);
            push @records, {
                direction => $dir,
                strain => $strain,
                tag => $tag,
                document => $actual_seed . ".xsd",
                bands => $bands_file,
                total_energy_kcalmol => $total_energy_kcalmol,
                total_energy_eV => $total_energy_eV,
                a => $a,
                b => $b,
                c => $c,
                v => $v,
                v0 => $v0
            };
        }
    }

    return @records;
}

sub make_strained_copy {
    my ($source, $name, $dir, $strain) = @_;

    my $doc = $source->SaveAs($name);
    $doc->MakeP1 if $MakeP1BeforeStrain =~ /^y/i;

    my $atoms = $doc->AsymmetricUnit->Atoms;
    my @frac;
    foreach my $atom (@$atoms) {
        my $f = $atom->FractionalXYZ;
        push @frac, Point(X => $f->X, Y => $f->Y, Z => $f->Z);
    }

    my ($a, $b, $c) = expected_lengths_after_strain($dir, $strain);
    my $lat = $doc->Lattice3D;
    $lat->LengthA = $a;
    $lat->LengthB = $b;
    $lat->LengthC = $c;

    for (my $i = 0; $i < @$atoms; ++$i) {
        $atoms->[$i]->FractionalXYZ = $frac[$i];
    }

    $doc->Save;
    return $doc;
}

sub result_total_energy {
    my ($results) = @_;
    my $energy_kcalmol = "";
    eval { $energy_kcalmol = $results->TotalEnergy; };
    return ("", "") if $@ || !defined($energy_kcalmol) || $energy_kcalmol eq "";
    $energy_kcalmol += 0.0;
    return ($energy_kcalmol, $energy_kcalmol * $KCALMOL_TO_EV_PER_CELL);
}

# ================= 3D k-grid and CASTEP Files input =================

sub build_mass_tensor_kpoint_grid {
    my ($data, $vbm, $cbm) = @_;
    my (%seen, @points);
    add_edge_kgrid(\@points, \%seen, $data->{kfrac}[$vbm->{kindex}]);
    add_edge_kgrid(\@points, \%seen, $data->{kfrac}[$cbm->{kindex}]);
    return @points;
}

sub add_edge_kgrid {
    my ($points, $seen, $center) = @_;
    for (my $ia = -$MassGridHalfWidth; $ia <= $MassGridHalfWidth; ++$ia) {
        for (my $ib = -$MassGridHalfWidth; $ib <= $MassGridHalfWidth; ++$ib) {
            for (my $ic = -$MassGridHalfWidth; $ic <= $MassGridHalfWidth; ++$ic) {
                my @k = (
                    wrap_frac01($center->[0] + $ia * $MassGridStepFrac),
                    wrap_frac01($center->[1] + $ib * $MassGridStepFrac),
                    wrap_frac01($center->[2] + $ic * $MassGridStepFrac)
                );
                my $key = sprintf("%.10f %.10f %.10f", @k);
                next if $seen->{$key}++;
                push @$points, \@k;
            }
        }
    }
}

sub spectral_kpoint_list_block {
    my ($points) = @_;
    die "No k-points were generated for mass tensor calculation.\n" unless @$points;
    my $w = 1.0 / scalar(@$points);
    my $txt = "%BLOCK SPECTRAL_KPOINT_LIST\n";
    foreach my $k (@$points) {
        $txt .= sprintf("  %.10f  %.10f  %.10f  %.14f\n", $k->[0], $k->[1], $k->[2], $w);
    }
    $txt .= "%ENDBLOCK SPECTRAL_KPOINT_LIST\n";
    return $txt;
}

sub make_mass_tensor_cell_text {
    my ($doc, $kblock) = @_;
    my $lat = $doc->Lattice3D;
    my @a = @{point_to_vec($lat->VectorA)};
    my @b = @{point_to_vec($lat->VectorB)};
    my @c = @{point_to_vec($lat->VectorC)};

    my $txt = "";
    $txt .= "%BLOCK LATTICE_CART\n";
    $txt .= "ang\n";
    $txt .= sprintf("  %.12f  %.12f  %.12f\n", @a);
    $txt .= sprintf("  %.12f  %.12f  %.12f\n", @b);
    $txt .= sprintf("  %.12f  %.12f  %.12f\n", @c);
    $txt .= "%ENDBLOCK LATTICE_CART\n\n";

    $txt .= "%BLOCK POSITIONS_FRAC\n";
    foreach my $atom (@{$doc->AsymmetricUnit->Atoms}) {
        my $f = $atom->FractionalXYZ;
        $txt .= sprintf("  %-4s  %.12f  %.12f  %.12f\n", $atom->ElementSymbol, $f->X, $f->Y, $f->Z);
    }
    $txt .= "%ENDBLOCK POSITIONS_FRAC\n\n";
    $txt .= "KPOINT_MP_SPACING $MassSCFKPointMPSpacing\n\n";
    $txt .= $kblock;
    return $txt;
}

sub make_mass_tensor_param_text {
    my $txt = "";
    $txt .= "task                    : spectral\n";
    $txt .= "spectral_task           : bandstructure\n";
    $txt .= "xc_functional           : $MassXCFunctional\n";
    $txt .= "cut_off_energy          : $MassCutOffEnergy\n";
    $txt .= "spectral_nextra_bands   : $CASTEPBandStructureSettings{BandStructureNumExtraBands}\n";
    $txt .= "spectral_eigenvalue_tol : $CASTEPBandStructureSettings{BandStructureEnergyTolerance} eV\n";
    $txt .= "write_bands             : true\n";
    return $txt;
}

# ================= .bands parser for reference edges =================

sub read_castep_bands {
    my ($file) = @_;
    open my $fh, "<", $file or die "Could not open '$file': $!\n";

    my $nk = last_number(next_nonempty_line($fh));
    my $nspin = last_number(next_nonempty_line($fh));
    my $nelec = last_number(next_nonempty_line($fh));
    my $nband = last_number(next_nonempty_line($fh));
    my $fermi_ha = last_number(next_nonempty_line($fh));

    my $cell_header = next_nonempty_line($fh);
    die "Expected 'Unit cell vectors' line, got: $cell_header\n"
        unless $cell_header =~ /unit\s+cell\s+vectors/i;
    for (my $i = 0; $i < 3; ++$i) { next_nonempty_line($fh); }

    my (@kfrac, @eig_eV);
    for (my $read = 0; $read < $nk; ++$read) {
        my $line = next_nonempty_line($fh);
        die "Expected K-point line, got: $line\n" unless $line =~ /k-?point/i;
        my @tok = split(/\s+/, trim($line));
        my $ik = int($tok[1]) - 1;
        $kfrac[$ik] = [$tok[2] + 0.0, $tok[3] + 0.0, $tok[4] + 0.0];

        for (my $is = 0; $is < $nspin; ++$is) {
            my $spin_line = next_nonempty_line($fh);
            die "Expected Spin component line, got: $spin_line\n" unless $spin_line =~ /spin\s+component/i;
            for (my $ib = 0; $ib < $nband; ++$ib) {
                my @nums = numbers_from_line(next_nonempty_line($fh));
                die "Could not read eigenvalue.\n" unless @nums;
                $eig_eV[$is][$ik][$ib] = $nums[-1] * $HARTREE_TO_EV;
            }
        }
    }
    close $fh;

    return {
        nk => int($nk),
        nspin => int($nspin),
        nband => int($nband),
        fermi_eV => $fermi_ha * $HARTREE_TO_EV,
        kfrac => \@kfrac,
        eig => \@eig_eV
    };
}

sub find_band_edges {
    my ($data) = @_;
    my ($vbm, $cbm);
    my $ef = $data->{fermi_eV};
    for (my $is = 0; $is < $data->{nspin}; ++$is) {
        for (my $ik = 0; $ik < $data->{nk}; ++$ik) {
            for (my $ib = 0; $ib < $data->{nband}; ++$ib) {
                my $e = $data->{eig}[$is][$ik][$ib];
                if ($e <= $ef + $FermiTolerance_eV) {
                    $vbm = { spin => $is, kindex => $ik, band => $ib, energy_eV => $e }
                        if !defined($vbm) || $e > $vbm->{energy_eV};
                }
                if ($e >= $ef - $FermiTolerance_eV) {
                    $cbm = { spin => $is, kindex => $ik, band => $ib, energy_eV => $e }
                        if !defined($cbm) || $e < $cbm->{energy_eV};
                }
            }
        }
    }
    return ($vbm, $cbm);
}

sub next_nonempty_line {
    my ($fh) = @_;
    while (defined(my $line = <$fh>)) {
        $line =~ s/\r?\n$//;
        next if $line =~ /^\s*$/;
        return $line;
    }
    die "Unexpected end of file.\n";
}

sub numbers_from_line {
    my ($line) = @_;
    my @nums = ($line =~ /[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?/g);
    foreach my $v (@nums) { $v =~ tr/Dd/Ee/; }
    return @nums;
}

sub last_number {
    my ($line) = @_;
    my @nums = numbers_from_line($line);
    die "No number found in line: $line\n" unless @nums;
    return $nums[-1];
}

sub trim {
    my ($s) = @_;
    $s =~ s/^\s+//;
    $s =~ s/\s+$//;
    return $s;
}

# ================= Geometry and file helpers =================

sub expected_cell_after_strain {
    my ($dir, $strain) = @_;
    my ($a, $b, $c) = expected_lengths_after_strain($dir, $strain);
    return ($a, $b, $c, $v0 * (1.0 + $strain));
}

sub expected_lengths_after_strain {
    my ($dir, $strain) = @_;
    my ($a, $b, $c) = ($a0, $b0, $c0);
    if ($dir =~ /^A$/i) {
        $a *= 1.0 + $strain;
    } elsif ($dir =~ /^B$/i) {
        $b *= 1.0 + $strain;
    } elsif ($dir =~ /^C$/i) {
        $c *= 1.0 + $strain;
    } else {
        die "Unknown direction '$dir'. Use A, B, or C.\n";
    }
    return ($a, $b, $c);
}

sub find_band_file {
    my ($seed, $keyword, $dirs, $min_mtime) = @_;
    my @roots = expanded_search_dirs($dirs);
    my (%seen, @found);
    foreach my $dir (@roots) {
        next unless -d $dir;
        find({
            wanted => sub {
                return unless -f $_;
                return unless $_ =~ /\.bands$/i;
                return unless $_ =~ /\Q$keyword\E/i;
                my $path = normalize_path($File::Find::name);
                return if defined($min_mtime) && (stat($path))[9] < $min_mtime;
                return unless $seed eq "" || file_name($path) =~ /\Q$seed\E/i || dirname($path) =~ /\Q$seed\E/i;
                return if $seen{lc($path)}++;
                push @found, $path;
            },
            no_chdir => 1
        }, $dir);
    }
    return undef unless @found;
    @found = sort { (stat($b))[9] <=> (stat($a))[9] } @found;
    return $found[0];
}

sub find_project_document_folder {
    my ($seed, $dirs) = @_;
    my @roots = expanded_search_dirs($dirs);
    my (%seen, @found);
    foreach my $dir (@roots) {
        next unless -d $dir;
        find({
            wanted => sub {
                return unless -f $_;
                return unless $_ =~ /^\Q$seed\E\.xsd$/i;
                my $path = normalize_path($File::Find::name);
                return if $seen{lc($path)}++;
                push @found, dirname($path);
            },
            no_chdir => 1
        }, $dir);
    }
    return undef unless @found;
    @found = sort { (stat($b))[9] <=> (stat($a))[9] } @found;
    return $found[0];
}

sub expanded_search_dirs {
    my ($dirs) = @_;
    my (%seen, @out);
    foreach my $dir (getcwd(), @$dirs) {
        $dir = normalize_path($dir);
        next if $dir eq "" || $seen{lc($dir)}++;
        push @out, $dir if -d $dir;
    }
    return @out;
}

sub write_text_document {
    my ($name, $content) = @_;
    my $path = normalize_path(getcwd() . "/" . $name);
    write_external_text_file($path, $content);
    eval { Documents->Import($path); };
}

sub write_external_text_file {
    my ($path, $content) = @_;
    open my $fh, ">", $path or die "Could not write '$path': $!\n";
    print $fh $content;
    close $fh;
}

sub point_to_vec {
    my ($p) = @_;
    return [$p->X, $p->Y, $p->Z];
}

sub wrap_frac01 {
    my ($x) = @_;
    $x -= int($x);
    $x += 1.0 if $x < 0.0;
    return $x;
}

sub strip_extension {
    my ($name) = @_;
    $name =~ s/\.[^.]+$//;
    return $name;
}

sub strain_tag {
    my ($strain) = @_;
    my $prefix = $strain < -1.0e-12 ? "m" : "p";
    my $milli = int(abs($strain) * 1000.0 + 0.5);
    return $prefix . sprintf("%03d", $milli);
}

sub short_seed {
    my ($seed) = @_;
    $seed =~ s/[^A-Za-z0-9]+//g;
    $seed = "Seed" if $seed eq "";
    return substr($seed, 0, 15);
}

sub short_deformation_seed {
    my ($base, $dir, $tag) = @_;
    $base =~ s/[^A-Za-z0-9]+//g;
    $base = "S" if $base eq "";
    my $seed = substr($base, 0, 4) . "DP" . $dir . $tag;
    return substr($seed, 0, 15);
}

sub optimized_seed_from_doc {
    my ($doc) = @_;
    my $seed = $doc->Name . $InitialOptSeedSuffix;
    $seed =~ s/[^A-Za-z0-9]+//g;
    $seed = "OPT" if $seed eq "";
    return substr($seed, 0, 15);
}

sub file_name {
    my ($path) = @_;
    $path = normalize_path($path);
    $path =~ s{.*/}{};
    return $path;
}

sub dirname {
    my ($path) = @_;
    $path = normalize_path($path);
    $path =~ s{/[^/]*$}{};
    return $path;
}

sub normalize_path {
    my ($path) = @_;
    return "" unless defined $path;
    $path =~ s/\\/\//g;
    $path =~ s/\/+/\//g;
    return $path;
}

# ================= Manifest output =================

sub write_linked_calc_std {
    my ($mass_info, $records) = @_;
    my $std = Documents->New($OutputStdName);
    my $r = 0;

    $std->Cell($r, 0) = "InputDocument";       $std->Cell($r++, 1) = $InputDocumentName;
    $std->Cell($r, 0) = "WorkingDocument";     $std->Cell($r++, 1) = $source_doc->Name . ".xsd";
    $std->Cell($r, 0) = "Temperature_K";       $std->Cell($r++, 1) = $TemperatureK;
    $std->Cell($r, 0) = "ReferenceBandsFile";  $std->Cell($r++, 1) = $mass_info->{reference_bands};
    $std->Cell($r, 0) = "MassTensorBandsFile"; $std->Cell($r++, 1) = $mass_info->{mass_bands};
    $std->Cell($r, 0) = "MassTensorSeed";      $std->Cell($r++, 1) = $mass_info->{mass_seed};
    $std->Cell($r, 0) = "MassGridHalfWidth";   $std->Cell($r++, 1) = $MassGridHalfWidth;
    $std->Cell($r, 0) = "MassGridStepFrac";    $std->Cell($r++, 1) = $MassGridStepFrac;
    $std->Cell($r, 0) = "GeneratedKPoints";    $std->Cell($r++, 1) = $mass_info->{generated_kpoints};
    $std->Cell($r, 0) = "V0_A3";               $std->Cell($r++, 1) = $v0;
    $std->Cell($r, 0) = "Reference_VBM_eV";    $std->Cell($r++, 1) = $mass_info->{vbm_eV};
    $std->Cell($r, 0) = "Reference_CBM_eV";    $std->Cell($r++, 1) = $mass_info->{cbm_eV};

    $r++;
    $std->Cell($r++, 0) = "DeformationRecords";
    my @header = qw(Direction Strain Tag Document BandsFile TotalEnergy_kcalmol TotalEnergy_eV a_A b_A c_A V_A3 V0_A3);
    for (my $c = 0; $c < @header; ++$c) { $std->Cell($r, $c) = $header[$c]; }
    ++$r;

    foreach my $rec (@$records) {
        $std->Cell($r, 0)  = $rec->{direction};
        $std->Cell($r, 1)  = $rec->{strain};
        $std->Cell($r, 2)  = $rec->{tag};
        $std->Cell($r, 3)  = $rec->{document};
        $std->Cell($r, 4)  = $rec->{bands};
        $std->Cell($r, 5)  = $rec->{total_energy_kcalmol};
        $std->Cell($r, 6)  = $rec->{total_energy_eV};
        $std->Cell($r, 7)  = $rec->{a};
        $std->Cell($r, 8)  = $rec->{b};
        $std->Cell($r, 9)  = $rec->{c};
        $std->Cell($r, 10) = $rec->{v};
        $std->Cell($r, 11) = $rec->{v0};
        ++$r;
    }
}
