#!perl
use strict;
use warnings;
use MaterialsScript qw(:all);

# -------------------------------------------------------------------------
# Forcite Fe plate heating + water evaporation model.
#
# Model assumptions:
#   - Input is an .xsd containing an Fe slab and water molecules.
#   - Fe slab has about 10 atomic layers along SurfaceAxis.
#   - Bottom layers 1-3 are treated as a fully insulated/support region:
#       fixed XYZ, zero velocity.
#   - Remaining Fe layers are the sustained hot plate.
#   - Water atoms are free. Evaporation is counted by water O atoms whose
#       coordinate is above the Fe surface by EvaporationHeight_A.
#
# Important physics note:
#   Forcite is classical MD. To keep the iron plate shape while maintaining a
#   hot surface, this script uses a controlled hot-substrate model:
#     1. Hot Fe layers are re-assigned Maxwell velocities at HotFeTemperature
#        before every short dynamics chunk.
#     2. If RestoreFeShapeEachChunk is Yes, Fe coordinates are restored to the
#        initial slab coordinates before each chunk.
#   This is a heat-bath boundary model, not a fully free high-temperature
#   metal deformation simulation.
# -------------------------------------------------------------------------

# ================= User input =================

# Empty means use active document. Otherwise use exact project document name.
my $InputDocumentName = "DF.xsd";

my $MakeWorkingCopy = "Yes";
my $WorkingDocumentName = "";
my $RunId = int(rand(9000)) + 1000;

# Geometry and atom selection.
my $FeElement = "Fe";
my $WaterOElement = "O";
my $WaterHElement = "H";
my $SurfaceAxis = "Z";
my $InsulatedFrom = "Bottom";
my $InsulatedLayerCount = 3;
my $ExpectedFeLayerCount = 10;       # Only used for diagnostics.
my $LayerTolerance_A = 0.45;
my $EvaporationFrom = "Top";
my $EvaporationHeight_A = 8.0;       # Water O above surface + this distance.

# Thermal model.
my $WaterInitialTemperature = 300.0;
my $HotFeTemperature = 1800.0;
my $RandomSeed = 24681357;
my $RestoreFeShapeEachChunk = "Yes";
my $RemoveWaterDrift = "Yes";

# Forcite settings. If typing fails, try Universal or COMPASSII.
my $Forcefield = "COMPASSIII";
my $Quality = "Medium";
my $ChargeAssignment = "Forcefield assigned";
my $AssignForcefieldTypes = "Yes";
my $WriteLevel = "Summary";

# Optional pre-minimization with insulated layers fixed.
my $RunGeometryOptimization = "No";
my $GeometryOptimizationSteps = 500;

# Dynamics.
my $ProductionSteps = 50000;
my $ChunkSteps = 250;                # Short chunks help keep Fe shape stable.
my $ProductionEnsemble = "NVE";      # Heating is imposed by hot Fe velocities.
my $TimeStep_fs = 1;
my $TrajectoryFrequency = 50;
my $DataFrequency = 50;
my $EnergyDeviation = 1.0e8;

# Outputs.
# Yes = save one continuous Forcite restart/append production trajectory.
my $CreateCombinedTrajectory = "Yes";
my $CombinedTrajectoryName = "heat-all.xtd";
my $AnalysisStdName = "";
my $CreateSets = "Yes";
my $SaveAllAtEnd = "Yes";
my $AnalysisStartRow = 10;

# Forcite Dynamics data-file output. These names follow the
# MaterialsScript API Forcite Dynamics Output Settings table.
my $WriteTemperatureData = "Yes";
my $WriteEnergyData = "Yes";
my $WritePressureData = "Yes";
my $WriteVolumeData = "Yes";
my $WriteDensityData = "Yes";
my $WriteStressData = "Yes";
my $WriteLatticeVectorData = "Yes";
my $WriteEnergyComponentData = "Yes";

Documents->SaveAllAtEnd = $SaveAllAtEnd;
srand($RandomSeed);

# ================= Main =================

die "InsulatedLayerCount must be >= 1.\n" unless $InsulatedLayerCount >= 1;
die "ProductionSteps must be >= 1.\n" unless $ProductionSteps >= 1;
die "ChunkSteps must be >= 1.\n" unless $ChunkSteps >= 1;
die "TimeStep_fs must be > 0.\n" unless $TimeStep_fs > 0.0;
die "HotFeTemperature must be > 0.\n" unless $HotFeTemperature > 0.0;
die "WaterInitialTemperature must be > 0.\n" unless $WaterInitialTemperature > 0.0;

my $source_doc = get_input_document($InputDocumentName);
validate_atomistic_document($source_doc, "input document");
my $doc = prepare_working_document($source_doc);
validate_atomistic_document($doc, "working document");
my $reference_positions = capture_positions($doc);

my $setup = classify_structure($doc, $CreateSets);
print_setup_summary($doc, $setup);

configure_forcite_common();
restore_and_control_fe($doc, $setup, $reference_positions);
assign_water_temperature($doc, $setup, $WaterInitialTemperature);

if ($RunGeometryOptimization =~ /^y/i) {
    print "\nRunning constrained Fe/water geometry optimization...\n";
    Modules->Forcite->GeometryOptimization->Run(
        $doc,
        Settings(
            CurrentForcefield => $Forcefield,
            Quality => $Quality,
            WriteLevel => $WriteLevel,
            AssignForcefieldTypes => $AssignForcefieldTypes,
            ChargeAssignment => $ChargeAssignment,
            MaxIterations => $GeometryOptimizationSteps,
            OptimizeCell => "No"
        )
    );
    $reference_positions = capture_positions($doc);
    $setup = classify_structure($doc, "No");
    restore_and_control_fe($doc, $setup, $reference_positions);
}

my $std = create_analysis_table($setup);
my $combined_doc = undef;

run_production($doc, $setup, $reference_positions, $std, \$combined_doc);

$std->Save;
eval { $combined_doc->Save; } if defined $combined_doc;
print "\nFinished Fe plate water evaporation simulation.\n";
print "Analysis table: " . $std->Name . "\n";
print "Combined trajectory: " . $combined_doc->Name . "\n" if defined $combined_doc;

# ================= Document setup =================

sub get_input_document {
    my ($name) = @_;

    if (defined($name) && $name ne "") {
        my $doc = eval { $Documents{$name}; };
        die "Could not find project document '$name'. Check exact XSD name.\n"
            unless defined $doc;
        return $doc;
    }

    my $doc = eval { Documents->ActiveDocument; };
    die "No InputDocumentName was set and no active document was available.\n"
        unless defined $doc;
    return $doc;
}

sub validate_atomistic_document {
    my ($doc, $label) = @_;

    die "Undefined $label.\n" unless defined $doc;

    my $ok = eval {
        my $atoms = $doc->UnitCell->Atoms;
        defined($atoms) && $atoms->Count >= 0;
    };

    die "The $label is not an atomistic XSD document. "
      . "Set \$InputDocumentName to the exact Fe/water .xsd name, "
      . "or make that .xsd the active document before running this script.\n"
        unless $ok;
}

sub prepare_working_document {
    my ($source) = @_;
    return $source unless $MakeWorkingCopy =~ /^y/i;

    my $name = $WorkingDocumentName;
    $name = short_name("FEW$RunId", "xsd") if !defined($name) || $name eq "";
    print "Creating working copy: $name\n";
    return $source->SaveAs($name);
}

sub configure_forcite_common {
    Modules->Forcite->ChangeSettings(
        Settings(
            CurrentForcefield => $Forcefield,
            Quality => $Quality,
            WriteLevel => $WriteLevel,
            AssignForcefieldTypes => $AssignForcefieldTypes,
            ChargeAssignment => $ChargeAssignment
        )
    );
}

# ================= Selection and controls =================

sub classify_structure {
    my ($doc, $create_sets) = @_;
    my @atoms = @{$doc->UnitCell->Atoms};
    die "No atoms were found in " . $doc->Name . ".\n" unless @atoms;

    my @fe_records;
    my @water_indices;
    my @water_o_indices;
    my @water_h_indices;

    for (my $i = 0; $i < @atoms; ++$i) {
        my $atom = $atoms[$i];
        my $element = $atom->ElementSymbol;

        if ($element eq $FeElement) {
            push @fe_records, { atom => $atom, index => $i };
            next;
        }

        push @water_indices, $i;
        push @water_o_indices, $i if $element eq $WaterOElement;
        push @water_h_indices, $i if $element eq $WaterHElement;
    }

    die "No Fe atoms were found. Check FeElement.\n" unless @fe_records;
    die "No water oxygen atoms were found. Check WaterOElement or model contents.\n"
        unless @water_o_indices;

    my @layers = build_layers(\@fe_records, $SurfaceAxis, $LayerTolerance_A);
    die "Only " . scalar(@layers) . " Fe layers were detected; cannot fix "
      . "$InsulatedLayerCount layer(s).\n"
        if @layers <= $InsulatedLayerCount;

    my @fixed_layers = select_side_layers(\@layers, $InsulatedFrom,
                                          $InsulatedLayerCount);
    my %fixed = map { $_->{index} => 1 } map { @{$_->{records}} } @fixed_layers;

    my @fixed_indices;
    my @hot_fe_indices;
    foreach my $record (@fe_records) {
        if ($fixed{$record->{index}}) {
            push @fixed_indices, $record->{index};
        } else {
            push @hot_fe_indices, $record->{index};
        }
    }

    die "Fixed Fe selection is empty.\n" unless @fixed_indices;
    die "Hot Fe selection is empty.\n" unless @hot_fe_indices;

    my $surface_coord = fe_surface_coordinate(\@layers);
    my $evap_threshold = evaporation_threshold($surface_coord);

    if ($create_sets =~ /^y/i) {
        create_named_set($doc, "FE_FIXED_L1_3", \@fixed_indices);
        create_named_set($doc, "FE_HOT_L4_10", \@hot_fe_indices);
        create_named_set($doc, "WATER_ALL", \@water_indices);
        create_named_set($doc, "WATER_O", \@water_o_indices);
    }

    return {
        atoms => \@atoms,
        fe_layers => \@layers,
        fixed_indices => \@fixed_indices,
        hot_fe_indices => \@hot_fe_indices,
        water_indices => \@water_indices,
        water_o_indices => \@water_o_indices,
        water_h_indices => \@water_h_indices,
        surface_coord => $surface_coord,
        evap_threshold => $evap_threshold
    };
}

sub build_layers {
    my ($records, $axis, $tol) = @_;
    my @sorted = sort {
        atom_coord($a->{atom}, $axis) <=> atom_coord($b->{atom}, $axis)
    } @$records;
    my @layers;

    foreach my $record (@sorted) {
        my $coord = atom_coord($record->{atom}, $axis);

        if (!@layers || abs($coord - $layers[-1]->{center}) > $tol) {
            push @layers, { center => $coord, records => [$record] };
            next;
        }

        push @{$layers[-1]->{records}}, $record;
        my $n = scalar @{$layers[-1]->{records}};
        $layers[-1]->{center} =
            (($layers[-1]->{center} * ($n - 1)) + $coord) / $n;
    }

    return @layers;
}

sub select_side_layers {
    my ($layers, $side, $count) = @_;
    my @selected;

    if ($side =~ /^bottom$/i) {
        for (my $i = 0; $i < $count; ++$i) {
            push @selected, $layers->[$i];
        }
        return @selected;
    }

    if ($side =~ /^top$/i) {
        for (my $i = @$layers - $count; $i < @$layers; ++$i) {
            push @selected, $layers->[$i];
        }
        return @selected;
    }

    die "InsulatedFrom must be Top or Bottom.\n";
}

sub fe_surface_coordinate {
    my ($layers) = @_;
    if ($EvaporationFrom =~ /^top$/i) {
        return $layers->[-1]->{center};
    }
    if ($EvaporationFrom =~ /^bottom$/i) {
        return $layers->[0]->{center};
    }
    die "EvaporationFrom must be Top or Bottom.\n";
}

sub evaporation_threshold {
    my ($surface_coord) = @_;
    if ($EvaporationFrom =~ /^top$/i) {
        return $surface_coord + $EvaporationHeight_A;
    }
    return $surface_coord - $EvaporationHeight_A;
}

sub capture_positions {
    my ($doc) = @_;
    my @atoms = @{$doc->UnitCell->Atoms};
    my @positions;
    for (my $i = 0; $i < @atoms; ++$i) {
        $positions[$i] = [$atoms[$i]->X, $atoms[$i]->Y, $atoms[$i]->Z];
    }
    return \@positions;
}

sub restore_and_control_fe {
    my ($doc, $setup, $reference) = @_;

    unfix_indices($doc, $setup->{fixed_indices});
    unfix_indices($doc, $setup->{hot_fe_indices});

    restore_positions($doc, $setup->{fixed_indices}, $reference);
    restore_positions($doc, $setup->{hot_fe_indices}, $reference)
        if $RestoreFeShapeEachChunk =~ /^y/i;

    my @fixed_atoms = atoms_by_indices($doc, $setup->{fixed_indices});
    foreach my $atom (@fixed_atoms) {
        $atom->Velocity = Point(X => 0.0, Y => 0.0, Z => 0.0);
        $atom->Fix("XYZ");
    }

    my @hot_atoms = atoms_by_indices($doc, $setup->{hot_fe_indices});
    foreach my $atom (@hot_atoms) {
        my ($vx, $vy, $vz) = random_velocity_vector($atom, $HotFeTemperature);
        $atom->Velocity = Point(X => $vx, Y => $vy, Z => $vz);
    }

    remove_center_of_mass_drift(\@hot_atoms);
}

sub unfix_indices {
    my ($doc, $indices) = @_;
    my @atoms = atoms_by_indices($doc, $indices);
    foreach my $atom (@atoms) {
        eval { $atom->Unfix("XYZ"); };
    }
}

sub restore_positions {
    my ($doc, $indices, $reference) = @_;
    my @atoms = @{$doc->UnitCell->Atoms};
    foreach my $idx (@$indices) {
        next unless defined $reference->[$idx];
        my ($x, $y, $z) = @{$reference->[$idx]};
        $atoms[$idx]->XYZ = Point(X => $x, Y => $y, Z => $z);
    }
}

sub assign_water_temperature {
    my ($doc, $setup, $temperature) = @_;
    my @water_atoms = atoms_by_indices($doc, $setup->{water_indices});
    foreach my $atom (@water_atoms) {
        my ($vx, $vy, $vz) = random_velocity_vector($atom, $temperature);
        $atom->Velocity = Point(X => $vx, Y => $vy, Z => $vz);
    }
    remove_center_of_mass_drift(\@water_atoms) if $RemoveWaterDrift =~ /^y/i;
}

# ================= Dynamics =================

sub dynamics_settings {
    my ($steps, $trajectory_restart, $append_trajectory) = @_;

    $trajectory_restart = "No" unless defined $trajectory_restart;
    $append_trajectory = "No" unless defined $append_trajectory;

    return Settings(
        CurrentForcefield => $Forcefield,
        Quality => $Quality,
        WriteLevel => $WriteLevel,
        AssignForcefieldTypes => $AssignForcefieldTypes,
        ChargeAssignment => $ChargeAssignment,
        Ensemble3D => $ProductionEnsemble,
        Ensemble2D => $ProductionEnsemble,
        Ensemble0D => $ProductionEnsemble,
        Temperature => $WaterInitialTemperature,
        InitialVelocities => "Current",
        RandomSeed => $RandomSeed,
        TimeStep => $TimeStep_fs,
        NumberOfSteps => $steps,
        TrajectoryFrequency => safe_frequency($steps, $TrajectoryFrequency),
        DataFrequency => safe_frequency($steps, $DataFrequency),
        WriteVelocities => "Yes",
        WriteForces => "Yes",
        WriteTemperatureData => $WriteTemperatureData,
        WriteEnergyData => $WriteEnergyData,
        WritePressureData => $WritePressureData,
        WriteVolumeData => $WriteVolumeData,
        WriteDensityData => $WriteDensityData,
        WriteStressData => $WriteStressData,
        WriteLatticeVectorData => $WriteLatticeVectorData,
        WriteEnergyComponentData => $WriteEnergyComponentData,
        EnergyDeviation => $EnergyDeviation,
        TrajectoryRestart => $trajectory_restart,
        AppendTrajectory => $append_trajectory
    );
}

sub run_production {
    my ($doc, $setup, $reference, $std, $combined_doc_ref) = @_;

    my $remaining = $ProductionSteps;
    my $cycle = 1;
    my $global_frame = 0;
    my $global_step_offset = 0;
    my $last_trajectory = undef;
    my $continuous_name = combined_trajectory_name();

    print "\nRunning Fe hot-plate / water evaporation production MD...\n";
    print "Total steps: $ProductionSteps; chunk steps: $ChunkSteps\n";
    print "Dynamics ensemble: $ProductionEnsemble; timestep: $TimeStep_fs fs\n";
    print "Continuous trajectory target: $continuous_name\n"
        if $CreateCombinedTrajectory =~ /^y/i;

    while ($remaining > 0) {
        set_current_frame_to_last($doc);

        my $steps_this = $remaining < $ChunkSteps ? $remaining : $ChunkSteps;
        my $frames_before = trajectory_frame_count($doc);

        $setup = classify_structure($doc, "No");
        restore_and_control_fe($doc, $setup, $reference);

        print "\nChunk $cycle: hot Fe = " . scalar(@{$setup->{hot_fe_indices}})
            . " atoms at $HotFeTemperature K; steps = $steps_this\n";

        my $trajectory_restart = defined($last_trajectory) ? "Yes" : "No";
        my $append_trajectory = defined($last_trajectory) ? "Yes" : "No";

        my $results = Modules->Forcite->Dynamics->Run(
            $doc,
            dynamics_settings($steps_this, $trajectory_restart, $append_trajectory)
        );
        my $result_summary = extract_forcite_result_summary($results);
        save_forcite_result_documents($results);

        my $traj_doc = eval { $results->Trajectory; };
        if (defined $traj_doc) {
            if ($cycle == 1 && $CreateCombinedTrajectory =~ /^y/i) {
                $traj_doc->Name = $continuous_name;
            }
            $last_trajectory = $traj_doc;
            $$combined_doc_ref = $traj_doc if $CreateCombinedTrajectory =~ /^y/i;
            print "Chunk trajectory: " . $traj_doc->Name . "\n";
            my $frames_after = trajectory_frame_count($traj_doc);
            my $first_new_frame = $frames_before + 1;
            my $last_new_frame = $frames_after;
            if ($first_new_frame <= $last_new_frame) {
                $global_frame = analyze_trajectory(
                    $traj_doc, $std, $setup, $cycle, $global_frame,
                    $global_step_offset,
                    safe_frequency($steps_this, $TrajectoryFrequency),
                    $result_summary,
                    $first_new_frame,
                    $last_new_frame
                );
            } else {
                print "Warning: no new trajectory frames were found for chunk $cycle.\n";
            }
        } else {
            print "Warning: Forcite returned no trajectory for chunk $cycle.\n";
        }

        $remaining -= $steps_this;
        $global_step_offset += $steps_this;
        last if $remaining <= 0;

        $doc = $last_trajectory if defined $last_trajectory;
        ++$cycle;
    }

    $doc = $last_trajectory if defined $last_trajectory;
    set_current_frame_to_last($doc);
    $setup = classify_structure($doc, "No");
    restore_and_control_fe($doc, $setup, $reference);
    eval { $doc->Save; };
    print "Final continuous trajectory frame Fe coordinates restored to the reference slab.\n";
}

sub combined_trajectory_name {
    my $name = $CombinedTrajectoryName;
    $name = "heat-all.xtd" if !defined($name) || $name eq "";
    $name = short_name("heatall$RunId", "xtd") if document_exists($name);
    return $name;
}

sub set_current_frame_to_last {
    my ($doc) = @_;

    my $ok = eval {
        my $trj = $doc->Trajectory;
        $trj->CurrentFrame = $trj->NumFrames if $trj->NumFrames > 0;
        1;
    };

    return $ok;
}

sub trajectory_frame_count {
    my ($doc) = @_;

    my $count = eval { $doc->Trajectory->NumFrames; };
    return 0 if $@ || !defined $count;
    return $count;
}

sub extract_forcite_result_summary {
    my ($results) = @_;

    return {
        Temperature => safe_result_scalar($results, "Temperature"),
        TotalEnergy => safe_result_scalar($results, "TotalEnergy"),
        PotentialEnergy => safe_result_scalar($results, "PotentialEnergy"),
        KineticEnergy => safe_result_scalar($results, "KineticEnergy"),
        Pressure => safe_result_scalar($results, "Pressure"),
        CellVolume => safe_result_scalar($results, "CellVolume"),
        Density => safe_result_scalar($results, "Density"),
        DataAsText => safe_result_document_name($results, "DataAsText"),
        DataAsChart => safe_result_document_name($results, "DataAsChart"),
        Report => safe_result_document_name($results, "Report"),
        EnergiesChart => safe_result_document_name($results, "EnergiesChart"),
        TemperatureChart => safe_result_document_name($results, "TemperatureChart")
    };
}

sub safe_result_scalar {
    my ($results, $property) = @_;
    return "" unless defined $results;

    my $value = eval { $results->$property; };
    return "" if $@ || !defined $value;
    return $value;
}

sub safe_result_document_name {
    my ($results, $property) = @_;
    return "" unless defined $results;

    my $doc = eval { $results->$property; };
    return "" if $@ || !defined $doc;

    my $name = eval { $doc->Name; };
    return "" if $@ || !defined $name;
    return $name;
}

sub save_forcite_result_documents {
    my ($results) = @_;
    return unless defined $results;

    foreach my $property (qw(
        Report EnergiesChart TemperatureChart DensityChart CellChart
        DataAsChart DataAsText ReducedTrajectory
    )) {
        my $doc = eval { $results->$property; };
        next if $@ || !defined $doc;
        eval { $doc->Save; };
    }
}

# ================= Analysis =================

sub create_analysis_table {
    my ($setup) = @_;
    my $name = $AnalysisStdName;
    $name = "FeWaterEvap_$RunId.std" if !defined($name) || $name eq "";

    my $std = Documents->New($name);
    my @headings = qw(
        Cycle LocalFrame GlobalFrame GlobalStep Time_ps
        EvaporatedWaterO TotalWaterO
        MaxOHeightAboveSurface_A AvgOHeightAboveSurface_A
        HotFeTemp_K WaterTemp_K SurfaceCoord_A EvapThreshold_A
        ForciteAvgTemperature_K ForciteAvgTotalEnergy_kcal_mol
        ForciteAvgPotentialEnergy_kcal_mol ForciteAvgKineticEnergy_kcal_mol
        ForciteAvgPressure_GPa ForciteAvgCellVolume_A3
        ForciteAvgDensity_g_cm3 ForciteDataAsText ForciteDataAsChart
        ForciteReport ForciteEnergiesChart ForciteTemperatureChart
    );

    for (my $c = 0; $c < @headings; ++$c) {
        $std->Cell(0, $c) = $headings[$c];
    }

    $std->Cell(2, 0) = "Setup";
    $std->Cell(3, 0) = "HotFeTemperature_K";
    $std->Cell(3, 1) = $HotFeTemperature;
    $std->Cell(4, 0) = "WaterInitialTemperature_K";
    $std->Cell(4, 1) = $WaterInitialTemperature;
    $std->Cell(5, 0) = "InsulatedLayerCount";
    $std->Cell(5, 1) = $InsulatedLayerCount;
    $std->Cell(6, 0) = "DetectedFeLayerCount";
    $std->Cell(6, 1) = scalar(@{$setup->{fe_layers}});
    $std->Cell(7, 0) = "EvaporationHeight_A";
    $std->Cell(7, 1) = $EvaporationHeight_A;
    $std->Cell(8, 0) = "DataFrequency_steps";
    $std->Cell(8, 1) = $DataFrequency;
    $std->Cell(9, 0) = "ForciteWriteDataFlags";
    $std->Cell(9, 1) = join(
        "; ",
        "Temperature=$WriteTemperatureData",
        "Energy=$WriteEnergyData",
        "Pressure=$WritePressureData",
        "Volume=$WriteVolumeData",
        "Density=$WriteDensityData",
        "Stress=$WriteStressData",
        "LatticeVector=$WriteLatticeVectorData",
        "EnergyComponent=$WriteEnergyComponentData"
    );

    return $std;
}

sub analyze_trajectory {
    my ($traj_doc, $std, $setup, $cycle, $global_frame,
        $global_step_offset, $traj_frequency, $result_summary,
        $start_frame, $end_frame) = @_;

    my $trajectory = eval { $traj_doc->Trajectory; };
    return $global_frame unless defined $trajectory;
    $result_summary = {} unless defined $result_summary;

    my $num_frames = $trajectory->NumFrames;
    $start_frame = 1 unless defined $start_frame;
    $end_frame = $num_frames unless defined $end_frame;
    $start_frame = 1 if $start_frame < 1;
    $end_frame = $num_frames if $end_frame > $num_frames;
    my $row = $AnalysisStartRow + $global_frame;

    for (my $frame = $start_frame; $frame <= $end_frame; ++$frame) {
        $trajectory->CurrentFrame = $frame;
        my $local_frame = $frame - $start_frame + 1;
        my $step = $global_step_offset + ($local_frame - 1) * $traj_frequency;
        my $time_ps = $step * $TimeStep_fs / 1000.0;

        my $evap = count_evaporated_water($traj_doc, $setup);
        my $hot_temp = kinetic_temperature(
            $traj_doc, $setup->{hot_fe_indices}
        );
        my $water_temp = kinetic_temperature(
            $traj_doc, $setup->{water_indices}
        );

        ++$global_frame;
        $std->Cell($row, 0) = $cycle;
        $std->Cell($row, 1) = $local_frame;
        $std->Cell($row, 2) = $global_frame;
        $std->Cell($row, 3) = $step;
        $std->Cell($row, 4) = $time_ps;
        $std->Cell($row, 5) = $evap->{count};
        $std->Cell($row, 6) = $evap->{total_o};
        $std->Cell($row, 7) = $evap->{max_height};
        $std->Cell($row, 8) = $evap->{avg_height};
        $std->Cell($row, 9) = $hot_temp;
        $std->Cell($row, 10) = $water_temp;
        $std->Cell($row, 11) = $setup->{surface_coord};
        $std->Cell($row, 12) = $setup->{evap_threshold};
        my @summary_keys = qw(
            Temperature TotalEnergy PotentialEnergy KineticEnergy
            Pressure CellVolume Density DataAsText DataAsChart
            Report EnergiesChart TemperatureChart
        );
        for (my $i = 0; $i < @summary_keys; ++$i) {
            my $value = $result_summary->{$summary_keys[$i]};
            $std->Cell($row, 13 + $i) = defined $value ? $value : "";
        }
        ++$row;
    }

    return $global_frame;
}

sub count_evaporated_water {
    my ($doc, $setup) = @_;
    my @atoms = @{$doc->UnitCell->Atoms};

    my $count = 0;
    my $total = 0;
    my $sum_height = 0.0;
    my $max_height = -1.0e30;

    foreach my $idx (@{$setup->{water_o_indices}}) {
        next if $idx < 0 || $idx >= @atoms;
        my $coord = atom_coord($atoms[$idx], $SurfaceAxis);
        my $height = ($EvaporationFrom =~ /^top$/i)
            ? $coord - $setup->{surface_coord}
            : $setup->{surface_coord} - $coord;

        ++$total;
        $sum_height += $height;
        $max_height = $height if $height > $max_height;

        if ($EvaporationFrom =~ /^top$/i) {
            ++$count if $coord >= $setup->{evap_threshold};
        } else {
            ++$count if $coord <= $setup->{evap_threshold};
        }
    }

    $max_height = 0.0 if $total == 0;
    my $avg_height = $total > 0 ? $sum_height / $total : 0.0;

    return {
        count => $count,
        total_o => $total,
        max_height => $max_height,
        avg_height => $avg_height
    };
}

sub kinetic_temperature {
    my ($doc, $indices) = @_;
    my @atoms = @{$doc->UnitCell->Atoms};
    my ($ke, $dof) = (0.0, 0);

    foreach my $idx (@$indices) {
        next if $idx < 0 || $idx >= @atoms;
        my $atom = $atoms[$idx];
        my $v = eval { $atom->Velocity; };
        next unless defined $v;
        my $m = $atom->Mass;
        next unless defined($m) && $m > 0.0;
        my $v2 = $v->X * $v->X + $v->Y * $v->Y + $v->Z * $v->Z;
        $ke += 0.5 * 0.00239005736153349 * $m * $v2;
        $dof += 3;
    }

    return 0.0 if $dof <= 0;
    return 2.0 * $ke / ($dof * 0.00198720425864083);
}

# ================= Utilities =================

sub print_setup_summary {
    my ($doc, $setup) = @_;
    my $layer_count = scalar(@{$setup->{fe_layers}});

    print "\nFe hot plate water evaporation setup\n";
    print "Working document: " . $doc->Name . "\n";
    print "Fe element: $FeElement\n";
    print "Detected Fe layers along $SurfaceAxis: $layer_count\n";
    print "Expected Fe layer count: $ExpectedFeLayerCount\n";
    print "Insulated fixed layers: $InsulatedLayerCount from $InsulatedFrom\n";
    print "Fixed Fe atoms: " . scalar(@{$setup->{fixed_indices}}) . "\n";
    print "Hot Fe atoms: " . scalar(@{$setup->{hot_fe_indices}}) . "\n";
    print "Water atoms: " . scalar(@{$setup->{water_indices}}) . "\n";
    print "Water O atoms: " . scalar(@{$setup->{water_o_indices}}) . "\n";
    print "Fe surface coordinate: " . $setup->{surface_coord} . " A\n";
    print "Evaporation threshold: " . $setup->{evap_threshold} . " A\n";
    print "Forcefield: $Forcefield\n";
    print "Hot Fe temperature: $HotFeTemperature K\n";
    print "Water initial temperature: $WaterInitialTemperature K\n";

    if ($ExpectedFeLayerCount > 0 && $layer_count != $ExpectedFeLayerCount) {
        print "Warning: detected layer count differs from ExpectedFeLayerCount. ";
        print "Adjust LayerTolerance_A if this is unexpected.\n";
    }
}

sub create_named_set {
    my ($doc, $name, $indices) = @_;
    return unless @$indices;

    delete_named_set($doc, $name);
    my @atoms = atoms_by_indices($doc, $indices);
    my $ok = eval {
        $doc->CreateSet($name, \@atoms);
        1;
    };
    print "Warning: could not create set '$name': $@\n" unless $ok;
}

sub delete_named_set {
    my ($doc, $name) = @_;
    my $set = eval { $doc->Sets($name); };
    eval { $set->Delete; } if defined $set;
}

sub atoms_by_indices {
    my ($doc, $indices) = @_;
    my @atoms = @{$doc->UnitCell->Atoms};
    my @selected;
    foreach my $idx (@$indices) {
        die "Atom index $idx is out of range in " . $doc->Name . ".\n"
            if $idx < 0 || $idx >= @atoms;
        push @selected, $atoms[$idx];
    }
    return @selected;
}

sub atom_coord {
    my ($atom, $axis) = @_;
    return $atom->X if $axis =~ /^x$/i;
    return $atom->Y if $axis =~ /^y$/i;
    return $atom->Z if $axis =~ /^z$/i;
    die "SurfaceAxis must be X, Y, or Z.\n";
}

sub random_velocity_vector {
    my ($atom, $temperature) = @_;

    my $mass = $atom->Mass;
    die "Atom mass is zero or undefined for " . $atom->ElementSymbol . ".\n"
        unless defined($mass) && $mass > 0.0;

    my $kb = 0.00198720425864083;
    my $conv = 0.00239005736153349;
    my $sigma = sqrt($kb * $temperature / ($conv * $mass));

    return (
        gaussian_random() * $sigma,
        gaussian_random() * $sigma,
        gaussian_random() * $sigma
    );
}

sub gaussian_random {
    my $u1 = rand();
    my $u2 = rand();
    $u1 = 1.0e-12 if $u1 < 1.0e-12;
    return sqrt(-2.0 * log($u1)) * cos(2.0 * 3.141592653589793 * $u2);
}

sub remove_center_of_mass_drift {
    my ($atoms) = @_;
    return unless @$atoms;

    my ($m_total, $px, $py, $pz) = (0.0, 0.0, 0.0, 0.0);
    foreach my $atom (@$atoms) {
        my $v = eval { $atom->Velocity; };
        next unless defined $v;
        my $m = $atom->Mass;
        next unless defined($m) && $m > 0.0;
        $m_total += $m;
        $px += $m * $v->X;
        $py += $m * $v->Y;
        $pz += $m * $v->Z;
    }
    return if $m_total <= 0.0;

    my ($dvx, $dvy, $dvz) = ($px / $m_total, $py / $m_total, $pz / $m_total);
    foreach my $atom (@$atoms) {
        my $v = eval { $atom->Velocity; };
        next unless defined $v;
        $atom->Velocity = Point(
            X => $v->X - $dvx,
            Y => $v->Y - $dvy,
            Z => $v->Z - $dvz
        );
    }
}

sub safe_frequency {
    my ($steps, $requested) = @_;
    return $steps if !defined($requested) || $requested < 1 || $requested > $steps;
    return $requested;
}

sub short_name {
    my ($prefix, $extension) = @_;
    $prefix =~ s/[^A-Za-z0-9]+//g;
    $prefix = "D" if $prefix eq "";
    my $max_base = 14 - length($extension);
    $prefix = substr($prefix, 0, $max_base) if length($prefix) > $max_base;

    for (my $attempt = 0; $attempt < 1000; ++$attempt) {
        my $suffix = $attempt == 0 ? "" : $attempt;
        my $base = substr($prefix . $suffix, 0, $max_base);
        my $name = "$base.$extension";
        next if document_exists($name);
        return $name;
    }

    die "Could not create a unique .$extension document name.\n";
}

sub document_exists {
    my ($name) = @_;
    my $doc = eval { $Documents{$name}; };
    return defined $doc;
}
