#!perl
use strict;
use warnings;
use Cwd qw(getcwd);
use File::Find;
use MaterialsScript qw(:all);

# -------------------------------------------------------------------------
# CASTEP mobility input workflow for a 3D bulk semiconductor.
# Author Jingsai Cheng
# E-mail jingsaicheng@gmail.com
# From Documents{"Si.xsd"}, this script calculates:
#   1. Directional deformation potentials along A, B, C:
#        E1_CBM_i = dE_CBM / d(epsilon_i)
#        E1_VBM_i = dE_VBM / d(epsilon_i)
#   2. Directional elastic moduli along A, B, C:
#        C_i = (1/V0) * d2Etotal / d(epsilon_i)^2
#   3. Electron and hole effective-mass tensors from an automatically generated
#      3D k-grid CASTEP .bands file around CBM/VBM.
#   4. Room-temperature mobility along A, B, C when all required data exist.
#
# Important:
#   Normal BandStru high-symmetry paths are not enough for a 3D effective-mass
#   tensor.  The script therefore generates SPECTRAL_KPOINT_LIST automatically.
# -------------------------------------------------------------------------

# ================= User input =================

my $InputDocumentName = "Si.xsd";

my @Directions = ("A", "B", "C");
my @Strains = (-0.010, -0.005, 0.000, 0.005, 0.010);

my $RunCASTEP = "Yes";
my $BandsNameKeyword = "BandStr";

# Restored previous local Materials Studio search paths.
my @BandsSearchDirs = (
    ".",
    getcwd(),
    "D:/msproject/bb_Files/Documents/test",

);

# Optional manual override for deformation-potential runs.
# Key format: "A:m010", "A:p000", "B:p005", ...
my %BandsFileByDirectionTag = (
    # "A:m010" => "D:/msproject/bb_Files/Documents/SiDPA_m010 CASTEP Energy/SiDPA_m010_BandStru.bands",
);

# Effective-mass tensor input.  If this is empty and $AutoMassTensorKGrid is
# Yes, the script first runs an unstrained BandStru calculation, locates VBM
# and CBM, generates a 3D SPECTRAL_KPOINT_LIST around them, runs CASTEP Files,
# and reads the resulting .bands.
my $MassTensorBandsFile = "";
my @MassTensorSearchDirs = @BandsSearchDirs;
my $AutoMassTensorKGrid = "Yes";
my $RunMassTensorFilesTask = "Yes";
my $MassTensorSeedPrefix = "SiM3D";
my $MassGridHalfWidth = 2;       # 2 gives 5x5x5 around each edge
my $MassGridStepFrac = 0.010;    # fractional reciprocal coordinate step
my $MassSCFKPointMPSpacing = "0.04 1/Ang";
my $MassCutOffEnergy = "500 eV";
my $MassXCFunctional = "PBE";

my $TemperatureK = 300.0;
my $OutputStdName = "ABC_Deformation_Elastic_Mass_Mobility.std";

my $MakeP1BeforeStrain = "Yes";
my $SaveAllAtEnd = "Yes";

# CASTEP settings for each strained structure.
my $CASTEPQuality = "Fine";
my $Pseudopotentials = "On the fly";
my $PropertiesKPointQuality = "Fine";
my $BandStructureEmax = 10;
my $BandStructureNumExtraBands = 12;
my $BandStructureEnergyTolerance = 1.0e-5;

# .bands parser settings.
my $CellVectorUnit = "bohr";       # CASTEP .bands cell vectors are normally bohr
my $FermiTolerance_eV = 1.0e-6;

# Effective-mass tensor fit settings.
my $MassTensorFitRadius_Ainv = 0.08;
my $MassTensorMaxEnergyWindow_eV = 1.0;
my $MassTensorMinimumFitPoints = 10;
my $UseFractionalMinimumImage = "Yes";

# ================= Constants =================

my $HARTREE_TO_EV = 27.211386245988;
my $BOHR_TO_A = 0.529177210903;
my $HBAR2_OVER_ME_EV_A2 = 7.619964231073;
my $KCALMOL_TO_EV_PER_CELL = 1.0 / 23.060548867;
my $EV_PER_A3_TO_GPA = 160.21766208;

my $PI = 3.14159265358979323846;
my $E_CHARGE = 1.602176634e-19;
my $KB = 1.380649e-23;
my $HBAR = 1.054571817e-34;
my $M0 = 9.1093837015e-31;

# ================= Main =================

Documents->SaveAllAtEnd = $SaveAllAtEnd;

my $source_doc = $Documents{$InputDocumentName};
die "Could not find Documents{\"$InputDocumentName\"}. Check the document name.\n"
    unless defined $source_doc;

my $lat0 = $source_doc->Lattice3D;
die "$InputDocumentName is not a 3D periodic document.\n" unless defined $lat0;

my $base_seed = strip_extension($InputDocumentName);
my $a0 = $lat0->LengthA;
my $b0 = $lat0->LengthB;
my $c0 = $lat0->LengthC;
my $v0 = $lat0->CellVolume;
my @real_dirs = (
    normalize_vec(point_to_vec($lat0->VectorA)),
    normalize_vec(point_to_vec($lat0->VectorB)),
    normalize_vec(point_to_vec($lat0->VectorC))
);

print "Source document: $InputDocumentName\n";
print "Directions: " . join(", ", @Directions) . "\n";
print "Temperature: $TemperatureK K\n";
print "Initial cell: a=$a0 A, b=$b0 A, c=$c0 A, V=$v0 A^3\n\n";

my $mass_result = calculate_mass_tensors_if_available();

my %dir_result;
foreach my $dir (@Directions) {
    print "================ Direction $dir ================\n";
    my @records = run_direction_series($dir);
    my $fit = fit_direction_results($dir, \@records);
    $dir_result{$dir} = $fit;
}

my @mobility_rows = calculate_mobility_rows(\%dir_result, $mass_result);

write_output_std(\%dir_result, $mass_result, \@mobility_rows);
print_summary(\%dir_result, $mass_result, \@mobility_rows);

# ================= Directional deformation/elastic series =================

sub run_direction_series {
    my ($dir) = @_;
    my @records;

    foreach my $strain (@Strains) {
        my $tag = strain_tag($strain);
        my $seed = short_seed($base_seed, $dir, $tag);
        my $xsd_name = "$seed.xsd";
        my $key = "$dir:$tag";

        print "--- $dir strain $strain ($tag) ---\n";

        my $work_doc = make_strained_copy($source_doc, $xsd_name, $dir, $strain);
        my $actual_seed = $work_doc->Name;
        print "Actual strained document: $actual_seed.xsd\n";

        my $bands_file;
        my $total_energy_kcalmol = "";
        my $total_energy_eV = "";

        if ($RunCASTEP =~ /^y/i) {
            print "Running CASTEP band structure for $actual_seed...\n";
            my $results = run_castep_band_structure($work_doc);
            ($total_energy_kcalmol, $total_energy_eV) = result_total_energy($results);

            $bands_file = find_band_file($actual_seed, $BandsNameKeyword, \@BandsSearchDirs);
            $bands_file = find_band_file($seed, $BandsNameKeyword, \@BandsSearchDirs)
                unless defined $bands_file;
            die "CASTEP finished, but no .bands file containing '$BandsNameKeyword' was found for $actual_seed.\n"
              . "Add the CASTEP project Documents folder to \@BandsSearchDirs.\n"
                unless defined $bands_file;
        } else {
            $bands_file = manual_deformation_bands_file($key);
            $bands_file = find_band_file($actual_seed, $BandsNameKeyword, \@BandsSearchDirs)
                unless defined $bands_file;
            $bands_file = find_band_file($seed, $BandsNameKeyword, \@BandsSearchDirs)
                unless defined $bands_file;
            die "No existing BandStru .bands found for $actual_seed and RunCASTEP is No.\n"
                unless defined $bands_file;
        }

        print "Using BandStru bands file: $bands_file\n";

        my $data = read_castep_bands($bands_file);
        my ($vbm, $cbm) = find_band_edges($data);
        die "Could not find VBM in $bands_file.\n" unless defined $vbm;
        die "Could not find CBM in $bands_file. The system may be metallic.\n" unless defined $cbm;

        my ($a, $b, $c, $v) = expected_cell_after_strain($dir, $strain);
        my $gap = $cbm->{energy_eV} - $vbm->{energy_eV};

        push @records, {
            direction => $dir,
            strain => $strain,
            tag => $tag,
            seed => $actual_seed,
            bands_file => $bands_file,
            total_energy_kcalmol => $total_energy_kcalmol,
            total_energy_eV => $total_energy_eV,
            a => $a,
            b => $b,
            c => $c,
            volume => $v,
            fermi_eV => $data->{fermi_eV},
            vbm => $vbm,
            cbm => $cbm,
            gap_eV => $gap
        };

        printf "VBM = %.8f eV, CBM = %.8f eV, gap = %.8f eV\n\n",
            $vbm->{energy_eV}, $cbm->{energy_eV}, $gap;
    }

    return @records;
}

sub fit_direction_results {
    my ($dir, $records) = @_;
    die "Direction $dir needs at least 3 strain points.\n" unless @$records >= 3;
    die "Direction $dir has missing total energies. Cannot fit elastic modulus.\n"
        unless all_records_have_total_energy($records);

    my @x = map { $_->{strain} } @$records;
    my @vbm = map { $_->{vbm}->{energy_eV} } @$records;
    my @cbm = map { $_->{cbm}->{energy_eV} } @$records;
    my @gap = map { $_->{gap_eV} } @$records;
    my @etot = map { $_->{total_energy_eV} } @$records;

    my $fit_vbm = linear_fit(\@x, \@vbm);
    my $fit_cbm = linear_fit(\@x, \@cbm);
    my $fit_gap = linear_fit(\@x, \@gap);
    my $fit_etot = quadratic_fit(\@x, \@etot);

    my $curvature_eV = $fit_etot->{curvature};
    my $C_GPa = $curvature_eV / $v0 * $EV_PER_A3_TO_GPA;

    printf "Direction %s: C = %.8f GPa, E1_CBM = %.8f eV, E1_VBM = %.8f eV\n\n",
        $dir, $C_GPa, $fit_cbm->{slope}, $fit_vbm->{slope};

    return {
        direction => $dir,
        records => $records,
        E1_VBM_fit => $fit_vbm,
        E1_CBM_fit => $fit_cbm,
        dEg_fit => $fit_gap,
        Etot_fit => $fit_etot,
        C_GPa => $C_GPa
    };
}

sub make_strained_copy {
    my ($source, $name, $dir, $strain) = @_;

    print "Creating strained copy: $name\n";
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

sub expected_cell_after_strain {
    my ($dir, $strain) = @_;
    my ($a, $b, $c) = expected_lengths_after_strain($dir, $strain);
    my $v = $v0;
    $v *= (1.0 + $strain) if $dir =~ /^[ABC]$/i;
    return ($a, $b, $c, $v);
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

sub run_castep_band_structure {
    my ($doc) = @_;
    my $results = Modules->CASTEP->Energy->Run($doc, Settings(
        Quality => $CASTEPQuality,
        Pseudopotentials => $Pseudopotentials,
        CalculateBandStructure => "Dispersion",
        PropertiesKPointQuality => $PropertiesKPointQuality,
        BandStructureEmax => $BandStructureEmax,
        BandStructureNumExtraBands => $BandStructureNumExtraBands,
        BandStructureEnergyTolerance => $BandStructureEnergyTolerance,
        BandStructureLineStyle => "Line"
    ));
    return $results;
}

sub result_total_energy {
    my ($results) = @_;
    my $energy_kcalmol = "";
    eval { $energy_kcalmol = $results->TotalEnergy; };
    return ("", "") if $@ || !defined($energy_kcalmol) || $energy_kcalmol eq "";
    $energy_kcalmol += 0.0;
    return ($energy_kcalmol, $energy_kcalmol * $KCALMOL_TO_EV_PER_CELL);
}

# ================= Effective-mass tensor =================

sub calculate_mass_tensors_if_available {
    my $file = resolve_mass_tensor_bands_file();
    if (!defined $file && $AutoMassTensorKGrid =~ /^y/i) {
        $file = auto_generate_and_run_mass_tensor_bands();
    }
    if (!defined $file) {
        print "Mass tensor .bands file was not provided/found. Tensor and mobility will be skipped.\n\n";
        return undef;
    }

    print "Reading effective-mass tensor .bands: $file\n";
    my $data = read_castep_bands($file);
    my @cell_A = cell_vectors_in_angstrom($data->{cell}, $CellVectorUnit);
    my @recip = reciprocal_vectors_2pi(\@cell_A);

    my ($vbm, $cbm) = find_band_edges($data);
    die "Could not find VBM in mass tensor .bands.\n" unless defined $vbm;
    die "Could not find CBM in mass tensor .bands.\n" unless defined $cbm;

    my $electron = fit_tensor_for_job($data, \@recip, make_job("electron", $cbm));
    my $hole = fit_tensor_for_job($data, \@recip, make_job("hole", $vbm));

    assign_directional_masses($electron);
    assign_directional_masses($hole);

    return {
        file => $file,
        data => $data,
        electron => $electron,
        hole => $hole
    };
}

sub auto_generate_and_run_mass_tensor_bands {
    print "No mass tensor .bands was provided. Building it automatically from unstrained structure.\n";

    my $ref_seed = short_ref_seed($base_seed);
    my $ref_doc = $source_doc->SaveAs("$ref_seed.xsd");
    my $actual_ref_seed = $ref_doc->Name;

    print "Running reference BandStru calculation for $actual_ref_seed...\n";
    run_castep_band_structure($ref_doc);
    my $ref_bands = find_band_file($actual_ref_seed, $BandsNameKeyword, \@BandsSearchDirs);
    $ref_bands = find_band_file($ref_seed, $BandsNameKeyword, \@BandsSearchDirs) unless defined $ref_bands;
    die "Reference BandStru .bands was not found after CASTEP run.\n"
        unless defined $ref_bands;

    print "Reference BandStru bands file: $ref_bands\n";
    my $ref_data = read_castep_bands($ref_bands);
    my ($vbm, $cbm) = find_band_edges($ref_data);
    die "Could not locate VBM/CBM from reference BandStru .bands.\n"
        unless defined $vbm && defined $cbm;

    printf "Reference VBM: energy %.8f eV, kpoint %d, band %d\n",
        $vbm->{energy_eV}, $vbm->{kindex} + 1, $vbm->{band} + 1;
    printf "Reference CBM: energy %.8f eV, kpoint %d, band %d\n",
        $cbm->{energy_eV}, $cbm->{kindex} + 1, $cbm->{band} + 1;

    my @kpoints = build_mass_tensor_kpoint_grid($ref_data, $vbm, $cbm);
    my $kblock = spectral_kpoint_list_block(\@kpoints);

    my $mass_seed = short_mass_seed($base_seed);
    my $mass_doc = $source_doc->SaveAs("$mass_seed.xsd");
    $mass_doc->MakeP1 if $MakeP1BeforeStrain =~ /^y/i;
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
        print "Also wrote CASTEP Files input directly to project folder: $project_dir\n";
    } else {
        print "Could not locate project folder for $actual_mass_seed.xsd before Files run.\n";
        print "If Files task ignores the generated .cell, manually place $actual_mass_seed.cell next to $actual_mass_seed.xsd.\n";
    }

    print "Generated CASTEP Files input: $actual_mass_seed.cell and $actual_mass_seed.param\n";
    print "Generated k-point block: ${actual_mass_seed}_SPECTRAL_KPOINT_LIST.txt\n";

    return undef unless $RunMassTensorFilesTask =~ /^y/i;

    print "Running CASTEP Files task for effective-mass tensor: $actual_mass_seed...\n";
    my $files_run_started = time;
    eval { Modules->CASTEP->Files->Run($mass_doc); };
    if ($@) {
        print "CASTEP Files task failed for $actual_mass_seed.\n";
        print "The input files were generated. Run $actual_mass_seed.xsd with CASTEP Files manually, then set MassTensorBandsFile to the resulting .bands.\n";
        print "Error was: $@\n";
        return undef;
    }

    my $mass_bands = find_band_file($actual_mass_seed, "", \@MassTensorSearchDirs, $files_run_started - 2);
    if (!defined $mass_bands) {
        print "CASTEP Files finished, but no new .bands file was found for $actual_mass_seed.\n";
        print "Set MassTensorBandsFile manually to the generated spectral .bands file.\n";
        return undef;
    }

    print "Mass tensor spectral .bands file: $mass_bands\n\n";
    return $mass_bands;
}

sub build_mass_tensor_kpoint_grid {
    my ($data, $vbm, $cbm) = @_;
    my %seen;
    my @points;
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

sub wrap_frac01 {
    my ($x) = @_;
    $x -= int($x);
    $x += 1.0 if $x < 0.0;
    return $x;
}

sub spectral_kpoint_list_block {
    my ($points) = @_;
    die "No k-points were generated for mass tensor calculation.\n" unless @$points;

    my $w = 1.0 / scalar(@$points);
    my $txt = "%BLOCK SPECTRAL_KPOINT_LIST\n";
    foreach my $k (@$points) {
        $txt .= sprintf("  %.10f  %.10f  %.10f  %.14f\n",
            $k->[0], $k->[1], $k->[2], $w);
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
        $txt .= sprintf("  %-4s  %.12f  %.12f  %.12f\n",
            $atom->ElementSymbol, $f->X, $f->Y, $f->Z);
    }
    $txt .= "%ENDBLOCK POSITIONS_FRAC\n\n";

    $txt .= "%BLOCK SPECIES_POT\n";
    $txt .= "%ENDBLOCK SPECIES_POT\n\n";
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
    $txt .= "spectral_nextra_bands   : $BandStructureNumExtraBands\n";
    $txt .= "spectral_eigenvalue_tol : $BandStructureEnergyTolerance eV\n";
    $txt .= "write_bands             : true\n";
    return $txt;
}

sub write_text_document {
    my ($name, $content) = @_;

    my $path = normalize_path(getcwd() . "/" . $name);
    open my $fh, ">", $path or die "Could not write '$path': $!\n";
    print $fh $content;
    close $fh;

    my $doc;
    eval { $doc = Documents->Import($path); };
    if ($@) {
        print "Wrote external file but could not import it into the MS project: $path\n";
        print "Import error: $@\n";
        return undef;
    }

    return $doc;
}

sub write_external_text_file {
    my ($path, $content) = @_;
    $path = normalize_path($path);
    open my $fh, ">", $path or die "Could not write '$path': $!\n";
    print $fh $content;
    close $fh;
}

sub fit_tensor_for_job {
    my ($data, $recip, $job) = @_;

    my $ik0 = $job->{kindex};
    my $k0_frac = $data->{kfrac}[$ik0];
    my $e0 = $job->{energy_eV};

    my (@rows, @energies, @used);
    for (my $ik = 0; $ik < $data->{nk}; ++$ik) {
        my $df = frac_delta($data->{kfrac}[$ik], $k0_frac);
        my $q = frac_to_cart($df, $recip);
        my $r = norm($q);
        next if $r > $MassTensorFitRadius_Ainv;

        my $e = $data->{eig}[$job->{spin}][$ik][$job->{band}];
        next if abs($e - $e0) > $MassTensorMaxEnergyWindow_eV;

        my ($qx, $qy, $qz) = @$q;
        push @rows, [
            1.0,
            $qx, $qy, $qz,
            0.5 * $qx * $qx,
            0.5 * $qy * $qy,
            0.5 * $qz * $qz,
            $qx * $qy,
            $qx * $qz,
            $qy * $qz
        ];
        push @energies, $e;
        push @used, $ik;
    }

    die "Not enough 3D k-points near $job->{carrier} edge: found "
      . scalar(@used) . ", need at least $MassTensorMinimumFitPoints.\n"
        if @used < $MassTensorMinimumFitPoints;

    my @coef = least_squares(\@rows, \@energies);
    die "Tensor fit failed for $job->{carrier}.\n" unless @coef == 10;

    my ($E0, $gx, $gy, $gz, $Hxx, $Hyy, $Hzz, $Hxy, $Hxz, $Hyz) = @coef;
    my @H = (
        [$Hxx, $Hxy, $Hxz],
        [$Hxy, $Hyy, $Hyz],
        [$Hxz, $Hyz, $Hzz]
    );

    my $carrier_sign = $job->{carrier} =~ /^hole/i ? -1.0 : 1.0;
    my @inv_mass = matrix_scale(\@H, $carrier_sign / $HBAR2_OVER_ME_EV_A2);
    my @mass_tensor = invert_3x3(\@inv_mass);
    my $r2 = fit_r2_matrix(\@rows, \@energies, \@coef);

    return {
        carrier => $job->{carrier},
        spin => $job->{spin},
        band => $job->{band},
        kindex => $job->{kindex},
        energy_eV => $job->{energy_eV},
        npoints => scalar(@used),
        used_kpoints => join(" ", map { $_ + 1 } @used),
        coeff => \@coef,
        H => \@H,
        inv_mass => \@inv_mass,
        mass_tensor => \@mass_tensor,
        r2 => $r2
    };
}

sub assign_directional_masses {
    my ($fit) = @_;
    my @labels = ("A", "B", "C");
    my %m;
    for (my $i = 0; $i < 3; ++$i) {
        my $u = $real_dirs[$i];
        my $inv = quadratic_form($fit->{inv_mass}, $u);
        $m{$labels[$i]} = abs($inv) > 1.0e-14 ? 1.0 / $inv : "";
    }
    $fit->{directional_mass} = \%m;
    if ($m{A} ne "" && $m{B} ne "" && $m{C} ne "") {
        $fit->{dos_mass} = (abs($m{A} * $m{B} * $m{C})) ** (1.0 / 3.0);
    } else {
        $fit->{dos_mass} = "";
    }
}

sub make_job {
    my ($carrier, $edge) = @_;
    return {
        carrier => $carrier,
        spin => $edge->{spin},
        band => $edge->{band},
        kindex => $edge->{kindex},
        energy_eV => $edge->{energy_eV}
    };
}

sub frac_delta {
    my ($kf, $k0) = @_;
    my @d = (
        $kf->[0] - $k0->[0],
        $kf->[1] - $k0->[1],
        $kf->[2] - $k0->[2]
    );
    if ($UseFractionalMinimumImage =~ /^y/i) {
        for (my $i = 0; $i < 3; ++$i) {
            $d[$i] -= nearest_integer($d[$i]);
        }
    }
    return \@d;
}

sub nearest_integer {
    my ($x) = @_;
    return int($x + 0.5) if $x >= 0;
    return -int(-$x + 0.5);
}

# ================= Mobility =================

sub calculate_mobility_rows {
    my ($dir_result, $mass_result) = @_;
    return () unless defined $mass_result;

    my @rows;
    foreach my $dir (@Directions) {
        my $d = $dir_result->{$dir};
        next unless defined $d;

        my $C_GPa = $d->{C_GPa};
        my $E1e = abs($d->{E1_CBM_fit}->{slope});
        my $E1h = abs($d->{E1_VBM_fit}->{slope});

        my $me_dir = abs($mass_result->{electron}->{directional_mass}->{$dir});
        my $mh_dir = abs($mass_result->{hole}->{directional_mass}->{$dir});
        my $me_dos = abs($mass_result->{electron}->{dos_mass});
        my $mh_dos = abs($mass_result->{hole}->{dos_mass});

        my $mu_e = mobility_3d_cm2Vs($C_GPa, $E1e, $me_dir, $me_dos);
        my $mu_h = mobility_3d_cm2Vs($C_GPa, $E1h, $mh_dir, $mh_dos);

        push @rows, {
            direction => $dir,
            C_GPa => $C_GPa,
            E1e_eV => $E1e,
            E1h_eV => $E1h,
            me_dir_m0 => $me_dir,
            mh_dir_m0 => $mh_dir,
            me_dos_m0 => $me_dos,
            mh_dos_m0 => $mh_dos,
            mu_e_cm2Vs => $mu_e,
            mu_h_cm2Vs => $mu_h
        };
    }
    return @rows;
}

sub mobility_3d_cm2Vs {
    my ($C_GPa, $E1_eV, $m_dir_m0, $m_dos_m0) = @_;
    return "" unless defined($C_GPa) && $C_GPa ne "";
    return "" unless defined($E1_eV) && abs($E1_eV) > 1.0e-12;
    return "" unless defined($m_dir_m0) && abs($m_dir_m0) > 1.0e-12;
    return "" unless defined($m_dos_m0) && abs($m_dos_m0) > 1.0e-12;

    my $C_Pa = $C_GPa * 1.0e9;
    my $E1_J = abs($E1_eV) * $E_CHARGE;
    my $m_dir = abs($m_dir_m0) * $M0;
    my $m_dos = abs($m_dos_m0) * $M0;

    my $mu_m2Vs =
        2.0 * sqrt(2.0 * $PI) * $E_CHARGE * ($HBAR ** 4) * $C_Pa
        / (3.0 * (($KB * $TemperatureK) ** 1.5) * $m_dir * ($m_dos ** 1.5) * ($E1_J ** 2));

    return $mu_m2Vs * 1.0e4;
}

# ================= .bands parser =================

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

    my @cell;
    for (my $i = 0; $i < 3; ++$i) {
        my @v = numbers_from_line(next_nonempty_line($fh));
        die "Could not read unit-cell vector.\n" unless @v >= 3;
        push @cell, [@v[0..2]];
    }

    my (@kfrac, @eig_eV);
    for (my $read = 0; $read < $nk; ++$read) {
        my $line = next_nonempty_line($fh);
        die "Expected K-point line, got: $line\n" unless $line =~ /k-?point/i;

        my @tok = split(/\s+/, trim($line));
        die "Malformed K-point line: $line\n" unless @tok >= 6;

        my $ik = int($tok[1]) - 1;
        $kfrac[$ik] = [$tok[2] + 0.0, $tok[3] + 0.0, $tok[4] + 0.0];

        for (my $is = 0; $is < $nspin; ++$is) {
            my $spin_line = next_nonempty_line($fh);
            die "Expected Spin component line, got: $spin_line\n"
                unless $spin_line =~ /spin\s+component/i;

            for (my $ib = 0; $ib < $nband; ++$ib) {
                my @nums = numbers_from_line(next_nonempty_line($fh));
                die "Could not read eigenvalue.\n" unless @nums;
                $eig_eV[$is][$ik][$ib] = $nums[-1] * $HARTREE_TO_EV;
            }
        }
    }

    close $fh;
    return {
        file => $file,
        nk => int($nk),
        nspin => int($nspin),
        nelec => $nelec + 0.0,
        nband => int($nband),
        fermi_eV => $fermi_ha * $HARTREE_TO_EV,
        cell => \@cell,
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
                    if (!defined($vbm) || $e > $vbm->{energy_eV}) {
                        $vbm = { spin => $is, kindex => $ik, band => $ib, energy_eV => $e };
                    }
                }
                if ($e >= $ef - $FermiTolerance_eV) {
                    if (!defined($cbm) || $e < $cbm->{energy_eV}) {
                        $cbm = { spin => $is, kindex => $ik, band => $ib, energy_eV => $e };
                    }
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

# ================= Fitting and linear algebra =================

sub linear_fit {
    my ($x, $y) = @_;
    my $n = scalar(@$x);
    die "Need at least 2 points for linear fit.\n" if $n < 2;

    my ($sx, $sy, $sxx, $sxy) = (0.0, 0.0, 0.0, 0.0);
    for (my $i = 0; $i < $n; ++$i) {
        $sx  += $x->[$i];
        $sy  += $y->[$i];
        $sxx += $x->[$i] * $x->[$i];
        $sxy += $x->[$i] * $y->[$i];
    }
    my $den = $n * $sxx - $sx * $sx;
    die "Linear fit failed because all x values are identical.\n" if abs($den) < 1.0e-30;

    my $slope = ($n * $sxy - $sx * $sy) / $den;
    my $intercept = ($sy - $slope * $sx) / $n;
    my $r2 = line_r2($x, $y, $intercept, $slope);
    return { slope => $slope, intercept => $intercept, r2 => $r2 };
}

sub quadratic_fit {
    my ($x, $y) = @_;
    my $n = scalar(@$x);
    die "Need at least 3 points for quadratic fit.\n" if $n < 3;

    my ($s0, $s1, $s2, $s3, $s4) = ($n, 0.0, 0.0, 0.0, 0.0);
    my ($t0, $t1, $t2) = (0.0, 0.0, 0.0);
    for (my $i = 0; $i < $n; ++$i) {
        my $xi = $x->[$i];
        my $yi = $y->[$i];
        my $x2 = $xi * $xi;
        $s1 += $xi;
        $s2 += $x2;
        $s3 += $x2 * $xi;
        $s4 += $x2 * $x2;
        $t0 += $yi;
        $t1 += $xi * $yi;
        $t2 += $x2 * $yi;
    }

    my @A = (
        [$s0, $s1, $s2, $t0],
        [$s1, $s2, $s3, $t1],
        [$s2, $s3, $s4, $t2]
    );
    my ($a, $b, $c) = solve_augmented(\@A);
    die "Quadratic fit failed.\n" unless defined $c;
    return {
        a => $a,
        b => $b,
        c => $c,
        curvature => 2.0 * $c,
        r2 => quadratic_r2($x, $y, $a, $b, $c)
    };
}

sub least_squares {
    my ($rows, $y) = @_;
    my $m = scalar(@{$rows->[0]});
    my @A;
    for (my $i = 0; $i < $m; ++$i) {
        for (my $j = 0; $j <= $m; ++$j) {
            $A[$i][$j] = 0.0;
        }
    }

    for (my $r = 0; $r < @$rows; ++$r) {
        for (my $i = 0; $i < $m; ++$i) {
            for (my $j = 0; $j < $m; ++$j) {
                $A[$i][$j] += $rows->[$r][$i] * $rows->[$r][$j];
            }
            $A[$i][$m] += $rows->[$r][$i] * $y->[$r];
        }
    }

    return solve_augmented(\@A);
}

sub solve_augmented {
    my ($A) = @_;
    my $n = scalar(@$A);
    for (my $col = 0; $col < $n; ++$col) {
        my $pivot = $col;
        for (my $r = $col + 1; $r < $n; ++$r) {
            $pivot = $r if abs($A->[$r][$col]) > abs($A->[$pivot][$col]);
        }
        return () if abs($A->[$pivot][$col]) < 1.0e-30;

        if ($pivot != $col) {
            my $tmp = $A->[$col];
            $A->[$col] = $A->[$pivot];
            $A->[$pivot] = $tmp;
        }

        my $p = $A->[$col][$col];
        for (my $j = $col; $j <= $n; ++$j) {
            $A->[$col][$j] /= $p;
        }

        for (my $r = 0; $r < $n; ++$r) {
            next if $r == $col;
            my $f = $A->[$r][$col];
            for (my $j = $col; $j <= $n; ++$j) {
                $A->[$r][$j] -= $f * $A->[$col][$j];
            }
        }
    }

    my @x;
    for (my $i = 0; $i < $n; ++$i) {
        push @x, $A->[$i][$n];
    }
    return @x;
}

sub invert_3x3 {
    my ($m) = @_;
    my $det =
        $m->[0][0] * ($m->[1][1] * $m->[2][2] - $m->[1][2] * $m->[2][1])
      - $m->[0][1] * ($m->[1][0] * $m->[2][2] - $m->[1][2] * $m->[2][0])
      + $m->[0][2] * ($m->[1][0] * $m->[2][1] - $m->[1][1] * $m->[2][0]);
    die "Cannot invert singular 3x3 tensor.\n" if abs($det) < 1.0e-30;

    return (
        [
            ($m->[1][1] * $m->[2][2] - $m->[1][2] * $m->[2][1]) / $det,
            ($m->[0][2] * $m->[2][1] - $m->[0][1] * $m->[2][2]) / $det,
            ($m->[0][1] * $m->[1][2] - $m->[0][2] * $m->[1][1]) / $det
        ],
        [
            ($m->[1][2] * $m->[2][0] - $m->[1][0] * $m->[2][2]) / $det,
            ($m->[0][0] * $m->[2][2] - $m->[0][2] * $m->[2][0]) / $det,
            ($m->[0][2] * $m->[1][0] - $m->[0][0] * $m->[1][2]) / $det
        ],
        [
            ($m->[1][0] * $m->[2][1] - $m->[1][1] * $m->[2][0]) / $det,
            ($m->[0][1] * $m->[2][0] - $m->[0][0] * $m->[2][1]) / $det,
            ($m->[0][0] * $m->[1][1] - $m->[0][1] * $m->[1][0]) / $det
        ]
    );
}

sub line_r2 {
    my ($x, $y, $a, $b) = @_;
    my $n = scalar(@$y);
    my $mean = 0.0;
    foreach my $yy (@$y) { $mean += $yy; }
    $mean /= $n;
    my ($ss_res, $ss_tot) = (0.0, 0.0);
    for (my $i = 0; $i < $n; ++$i) {
        my $pred = $a + $b * $x->[$i];
        $ss_res += ($y->[$i] - $pred) ** 2;
        $ss_tot += ($y->[$i] - $mean) ** 2;
    }
    return $ss_tot > 0.0 ? 1.0 - $ss_res / $ss_tot : 1.0;
}

sub quadratic_r2 {
    my ($x, $y, $a, $b, $c) = @_;
    my $n = scalar(@$y);
    my $mean = 0.0;
    foreach my $yy (@$y) { $mean += $yy; }
    $mean /= $n;
    my ($ss_res, $ss_tot) = (0.0, 0.0);
    for (my $i = 0; $i < $n; ++$i) {
        my $pred = $a + $b * $x->[$i] + $c * $x->[$i] * $x->[$i];
        $ss_res += ($y->[$i] - $pred) ** 2;
        $ss_tot += ($y->[$i] - $mean) ** 2;
    }
    return $ss_tot > 0.0 ? 1.0 - $ss_res / $ss_tot : 1.0;
}

sub fit_r2_matrix {
    my ($rows, $y, $coef) = @_;
    my $n = scalar(@$y);
    my $mean = 0.0;
    foreach my $yy (@$y) { $mean += $yy; }
    $mean /= $n;
    my ($ss_res, $ss_tot) = (0.0, 0.0);
    for (my $r = 0; $r < $n; ++$r) {
        my $pred = 0.0;
        for (my $i = 0; $i < @$coef; ++$i) {
            $pred += $rows->[$r][$i] * $coef->[$i];
        }
        $ss_res += ($y->[$r] - $pred) ** 2;
        $ss_tot += ($y->[$r] - $mean) ** 2;
    }
    return $ss_tot > 0.0 ? 1.0 - $ss_res / $ss_tot : 1.0;
}

# ================= Geometry helpers =================

sub cell_vectors_in_angstrom {
    my ($cell, $unit) = @_;
    my $scale = ($unit =~ /^bohr$/i || $unit =~ /^au$/i) ? $BOHR_TO_A : 1.0;
    my @out;
    foreach my $v (@$cell) {
        push @out, [$v->[0] * $scale, $v->[1] * $scale, $v->[2] * $scale];
    }
    return @out;
}

sub reciprocal_vectors_2pi {
    my ($cell) = @_;
    my ($a, $b, $c) = @$cell;
    my $vol = dot($a, cross($b, $c));
    die "Unit cell volume is zero.\n" if abs($vol) < 1.0e-20;
    my $factor = 2.0 * $PI / $vol;
    return (
        scale_vec(cross($b, $c), $factor),
        scale_vec(cross($c, $a), $factor),
        scale_vec(cross($a, $b), $factor)
    );
}

sub frac_to_cart {
    my ($kf, $basis) = @_;
    my ($b1, $b2, $b3) = @$basis;
    return [
        $kf->[0] * $b1->[0] + $kf->[1] * $b2->[0] + $kf->[2] * $b3->[0],
        $kf->[0] * $b1->[1] + $kf->[1] * $b2->[1] + $kf->[2] * $b3->[1],
        $kf->[0] * $b1->[2] + $kf->[1] * $b2->[2] + $kf->[2] * $b3->[2]
    ];
}

sub point_to_vec {
    my ($p) = @_;
    return [$p->X, $p->Y, $p->Z];
}

sub normalize_vec {
    my ($v) = @_;
    my $n = norm($v);
    die "Cannot normalize zero vector.\n" if $n < 1.0e-30;
    return [$v->[0] / $n, $v->[1] / $n, $v->[2] / $n];
}

sub scale_vec {
    my ($v, $s) = @_;
    return [$v->[0] * $s, $v->[1] * $s, $v->[2] * $s];
}

sub matrix_scale {
    my ($m, $s) = @_;
    my @out;
    for (my $i = 0; $i < 3; ++$i) {
        for (my $j = 0; $j < 3; ++$j) {
            $out[$i][$j] = $m->[$i][$j] * $s;
        }
    }
    return @out;
}

sub quadratic_form {
    my ($m, $u) = @_;
    my $sum = 0.0;
    for (my $i = 0; $i < 3; ++$i) {
        for (my $j = 0; $j < 3; ++$j) {
            $sum += $u->[$i] * $m->[$i][$j] * $u->[$j];
        }
    }
    return $sum;
}

sub cross {
    my ($u, $v) = @_;
    return [
        $u->[1] * $v->[2] - $u->[2] * $v->[1],
        $u->[2] * $v->[0] - $u->[0] * $v->[2],
        $u->[0] * $v->[1] - $u->[1] * $v->[0]
    ];
}

sub dot {
    my ($u, $v) = @_;
    return $u->[0] * $v->[0] + $u->[1] * $v->[1] + $u->[2] * $v->[2];
}

sub norm {
    my ($v) = @_;
    return sqrt(dot($v, $v));
}

# ================= File resolving =================

sub manual_deformation_bands_file {
    my ($key) = @_;
    return undef unless exists $BandsFileByDirectionTag{$key};
    my $file = normalize_path($BandsFileByDirectionTag{$key});
    die "Manual bands file for $key does not exist: $file\n" unless -f $file;
    die "Manual bands file for $key does not contain '$BandsNameKeyword': $file\n"
        unless file_name($file) =~ /\Q$BandsNameKeyword\E/i;
    return $file;
}

sub resolve_mass_tensor_bands_file {
    if (defined($MassTensorBandsFile) && $MassTensorBandsFile ne "") {
        my $file = normalize_path($MassTensorBandsFile);
        return $file if -f $file;
        die "MassTensorBandsFile does not exist: $file\n";
    }
    return undef;
}

sub find_band_file {
    my ($seed, $keyword, $dirs, $min_mtime) = @_;
    my @roots = expanded_search_dirs($dirs);
    my %seen;
    my @found;

    foreach my $dir (@roots) {
        next unless -d $dir;
        find({
            wanted => sub {
                return unless -f $_;
                return unless $_ =~ /\.bands$/i;
                return unless $_ =~ /\Q$keyword\E/i;
                my $path = normalize_path($File::Find::name);
                return if defined($min_mtime) && (stat($path))[9] < $min_mtime;
                if (defined($seed) && $seed ne "") {
                    return unless file_name($path) =~ /\Q$seed\E/i
                               || dirname($path) =~ /\Q$seed\E/i;
                }
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
    my %seen;
    my @found;

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

sub portable_search_dirs {
    my @extra = @_;
    my @base = (
        ".",
        getcwd(),
        "Documents",
        "../Documents",
        "..",
        "../..",
        @extra
    );

    my %seen;
    my @out;
    foreach my $dir (@base) {
        $dir = normalize_path($dir);
        next if $dir eq "" || $seen{lc($dir)}++;
        push @out, $dir if -d $dir;
    }
    return @out;
}

sub expanded_search_dirs {
    my ($dirs) = @_;
    my %seen;
    my @out;
    foreach my $dir (@$dirs) {
        $dir = normalize_path($dir);
        next if $dir eq "" || $seen{lc($dir)}++;
        push @out, $dir if -d $dir;
    }
    return @out;
}

sub file_name {
    my ($path) = @_;
    $path = normalize_path($path);
    $path =~ s{^.*/}{};
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
    my ($base, $dir, $tag) = @_;
    $base =~ s/[^A-Za-z0-9]+//g;
    $base = "S" if $base eq "";
    my $seed = substr($base, 0, 4) . "DP" . $dir . $tag;
    return substr($seed, 0, 15);
}

sub short_ref_seed {
    my ($base) = @_;
    $base =~ s/[^A-Za-z0-9]+//g;
    $base = "S" if $base eq "";
    my $seed = substr($base, 0, 4) . "REF";
    return substr($seed, 0, 15);
}

sub short_mass_seed {
    my ($base) = @_;
    my $seed = $MassTensorSeedPrefix;
    if (!defined($seed) || $seed eq "") {
        $base =~ s/[^A-Za-z0-9]+//g;
        $base = "S" if $base eq "";
        $seed = substr($base, 0, 4) . "M3D";
    }
    $seed =~ s/[^A-Za-z0-9]+//g;
    $seed = "M3D" if $seed eq "";
    return substr($seed, 0, 15);
}

sub all_records_have_total_energy {
    my ($records) = @_;
    return 0 unless @$records >= 3;
    foreach my $rec (@$records) {
        return 0 unless defined($rec->{total_energy_eV}) && $rec->{total_energy_eV} ne "";
    }
    return 1;
}

# ================= Output =================

sub write_output_std {
    my ($dir_result, $mass_result, $mobility_rows) = @_;
    my $std = Documents->New($OutputStdName);

    my $r = 0;
    $std->Cell($r, 0) = "InputDocument";   $std->Cell($r++, 1) = $InputDocumentName;
    $std->Cell($r, 0) = "Temperature_K";   $std->Cell($r++, 1) = $TemperatureK;
    $std->Cell($r, 0) = "Definition_E1";   $std->Cell($r++, 1) = "E1=dE_edge/dstrain";
    $std->Cell($r, 0) = "Definition_C";    $std->Cell($r++, 1) = "C=(1/V0)*d2Etot/dstrain2";
    $std->Cell($r, 0) = "MobilityModel";   $std->Cell($r++, 1) = "3D acoustic deformation potential";
    $r++;

    $std->Cell($r, 0) = "Direction";
    $std->Cell($r, 1) = "C_GPa";
    $std->Cell($r, 2) = "E1_CBM_eV";
    $std->Cell($r, 3) = "Abs_E1_CBM_eV";
    $std->Cell($r, 4) = "E1_VBM_eV";
    $std->Cell($r, 5) = "Abs_E1_VBM_eV";
    $std->Cell($r, 6) = "dEg_dstrain_eV";
    $std->Cell($r, 7) = "R2_C";
    $std->Cell($r, 8) = "R2_CBM";
    $std->Cell($r, 9) = "R2_VBM";
    $r++;

    foreach my $dir (@Directions) {
        my $d = $dir_result->{$dir};
        next unless defined $d;
        $std->Cell($r, 0) = $dir;
        $std->Cell($r, 1) = $d->{C_GPa};
        $std->Cell($r, 2) = $d->{E1_CBM_fit}->{slope};
        $std->Cell($r, 3) = abs($d->{E1_CBM_fit}->{slope});
        $std->Cell($r, 4) = $d->{E1_VBM_fit}->{slope};
        $std->Cell($r, 5) = abs($d->{E1_VBM_fit}->{slope});
        $std->Cell($r, 6) = $d->{dEg_fit}->{slope};
        $std->Cell($r, 7) = $d->{Etot_fit}->{r2};
        $std->Cell($r, 8) = $d->{E1_CBM_fit}->{r2};
        $std->Cell($r, 9) = $d->{E1_VBM_fit}->{r2};
        $r++;
    }

    $r += 2;
    $std->Cell($r, 0) = "Per-strain data";
    $r++;
    my @cols = (
        "Direction", "Strain", "Tag", "Seed",
        "a_A", "b_A", "c_A", "Volume_A3",
        "Etot_kcal_mol", "Etot_eV_cell",
        "Fermi_eV", "VBM_eV", "CBM_eV", "Gap_eV",
        "VBM_spin", "VBM_kpoint", "VBM_band",
        "CBM_spin", "CBM_kpoint", "CBM_band",
        "BandStru_bands_file"
    );
    for (my $c = 0; $c < @cols; ++$c) {
        $std->Cell($r, $c) = $cols[$c];
    }
    $r++;

    foreach my $dir (@Directions) {
        my $d = $dir_result->{$dir};
        next unless defined $d;
        foreach my $rec (@{$d->{records}}) {
            $std->Cell($r, 0) = $dir;
            $std->Cell($r, 1) = $rec->{strain};
            $std->Cell($r, 2) = $rec->{tag};
            $std->Cell($r, 3) = $rec->{seed};
            $std->Cell($r, 4) = $rec->{a};
            $std->Cell($r, 5) = $rec->{b};
            $std->Cell($r, 6) = $rec->{c};
            $std->Cell($r, 7) = $rec->{volume};
            $std->Cell($r, 8) = $rec->{total_energy_kcalmol};
            $std->Cell($r, 9) = $rec->{total_energy_eV};
            $std->Cell($r, 10) = $rec->{fermi_eV};
            $std->Cell($r, 11) = $rec->{vbm}->{energy_eV};
            $std->Cell($r, 12) = $rec->{cbm}->{energy_eV};
            $std->Cell($r, 13) = $rec->{gap_eV};
            $std->Cell($r, 14) = $rec->{vbm}->{spin} + 1;
            $std->Cell($r, 15) = $rec->{vbm}->{kindex} + 1;
            $std->Cell($r, 16) = $rec->{vbm}->{band} + 1;
            $std->Cell($r, 17) = $rec->{cbm}->{spin} + 1;
            $std->Cell($r, 18) = $rec->{cbm}->{kindex} + 1;
            $std->Cell($r, 19) = $rec->{cbm}->{band} + 1;
            $std->Cell($r, 20) = $rec->{bands_file};
            $r++;
        }
    }

    $r += 2;
    write_mass_output($std, \$r, $mass_result);
    $r += 2;
    write_mobility_output($std, \$r, $mobility_rows);
}

sub write_mass_output {
    my ($std, $rref, $mass_result) = @_;
    my $r = $$rref;
    $std->Cell($r, 0) = "Effective mass tensor";
    if (!defined $mass_result) {
        $std->Cell($r, 1) = "Skipped: set MassTensorBandsFile to a 3D k-grid .bands file.";
        $$rref = $r + 1;
        return;
    }

    $std->Cell($r++, 1) = $mass_result->{file};
    foreach my $carrier ("electron", "hole") {
        my $fit = $mass_result->{$carrier};
        $std->Cell($r, 0) = $carrier;
        $std->Cell($r, 1) = "edge_energy_eV";
        $std->Cell($r, 2) = $fit->{energy_eV};
        $std->Cell($r, 3) = "spin";
        $std->Cell($r, 4) = $fit->{spin} + 1;
        $std->Cell($r, 5) = "kpoint";
        $std->Cell($r, 6) = $fit->{kindex} + 1;
        $std->Cell($r, 7) = "band";
        $std->Cell($r, 8) = $fit->{band} + 1;
        $std->Cell($r, 9) = "npoints";
        $std->Cell($r, 10) = $fit->{npoints};
        $std->Cell($r, 11) = "R2";
        $std->Cell($r, 12) = $fit->{r2};
        $r++;

        $std->Cell($r, 0) = "${carrier}_directional_mass_m0";
        $std->Cell($r, 1) = "m_A";
        $std->Cell($r, 2) = $fit->{directional_mass}->{A};
        $std->Cell($r, 3) = "m_B";
        $std->Cell($r, 4) = $fit->{directional_mass}->{B};
        $std->Cell($r, 5) = "m_C";
        $std->Cell($r, 6) = $fit->{directional_mass}->{C};
        $std->Cell($r, 7) = "m_DOS";
        $std->Cell($r, 8) = $fit->{dos_mass};
        $r++;

        $std->Cell($r, 0) = "${carrier}_mass_tensor_m0";
        $r++;
        for (my $i = 0; $i < 3; ++$i) {
            for (my $j = 0; $j < 3; ++$j) {
                $std->Cell($r, $j) = $fit->{mass_tensor}->[$i][$j];
            }
            $r++;
        }

        $std->Cell($r, 0) = "${carrier}_inverse_mass_tensor_1_m0";
        $r++;
        for (my $i = 0; $i < 3; ++$i) {
            for (my $j = 0; $j < 3; ++$j) {
                $std->Cell($r, $j) = $fit->{inv_mass}->[$i][$j];
            }
            $r++;
        }
        $r++;
    }
    $$rref = $r;
}

sub write_mobility_output {
    my ($std, $rref, $mobility_rows) = @_;
    my $r = $$rref;
    $std->Cell($r++, 0) = "Room-temperature mobility";
    my @cols = (
        "Direction", "C_GPa",
        "E1e_CBM_abs_eV", "m_e_dir_m0", "m_e_DOS_m0", "mu_e_cm2_Vs",
        "E1h_VBM_abs_eV", "m_h_dir_m0", "m_h_DOS_m0", "mu_h_cm2_Vs"
    );
    for (my $c = 0; $c < @cols; ++$c) {
        $std->Cell($r, $c) = $cols[$c];
    }
    $r++;

    if (!@$mobility_rows) {
        $std->Cell($r, 0) = "Skipped: effective-mass tensor was not available.";
        $$rref = $r + 1;
        return;
    }

    foreach my $row (@$mobility_rows) {
        $std->Cell($r, 0) = $row->{direction};
        $std->Cell($r, 1) = $row->{C_GPa};
        $std->Cell($r, 2) = $row->{E1e_eV};
        $std->Cell($r, 3) = $row->{me_dir_m0};
        $std->Cell($r, 4) = $row->{me_dos_m0};
        $std->Cell($r, 5) = $row->{mu_e_cm2Vs};
        $std->Cell($r, 6) = $row->{E1h_eV};
        $std->Cell($r, 7) = $row->{mh_dir_m0};
        $std->Cell($r, 8) = $row->{mh_dos_m0};
        $std->Cell($r, 9) = $row->{mu_h_cm2Vs};
        $r++;
    }
    $$rref = $r;
}

sub print_summary {
    my ($dir_result, $mass_result, $mobility_rows) = @_;
    print "STD generated: $OutputStdName\n\n";
    foreach my $dir (@Directions) {
        my $d = $dir_result->{$dir};
        next unless defined $d;
        printf "%s: C=%.6f GPa, E1e=%.6f eV, E1h=%.6f eV\n",
            $dir, $d->{C_GPa}, $d->{E1_CBM_fit}->{slope}, $d->{E1_VBM_fit}->{slope};
    }

    if (defined $mass_result) {
        print "\nEffective mass tensor fitted from: $mass_result->{file}\n";
        foreach my $carrier ("electron", "hole") {
            my $m = $mass_result->{$carrier}->{directional_mass};
            printf "%s masses: mA=%.6f, mB=%.6f, mC=%.6f, mdos=%.6f\n",
                $carrier, $m->{A}, $m->{B}, $m->{C}, $mass_result->{$carrier}->{dos_mass};
        }
    }

    if (@$mobility_rows) {
        print "\nRoom-temperature mobility (cm2/V/s):\n";
        foreach my $row (@$mobility_rows) {
            printf "%s: electron=%.6f, hole=%.6f\n",
                $row->{direction}, $row->{mu_e_cm2Vs}, $row->{mu_h_cm2Vs};
        }
    }
}
