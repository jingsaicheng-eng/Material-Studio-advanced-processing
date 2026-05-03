#!perl
use strict;
use warnings;
use MaterialsScript qw(:all);

# -------------------------------------------------------------------------
# Local center-heating / heat-transfer MD for a Cu surface with Forcite.
#
# Model:
#   - Bottom 5 atomic layers are fixed.
#   - All other atoms start at 300 K.
#   - The heat source is a circular region at the surface center.
#   - Heat-source radius is 3 Angstrom.
#   - Heat-source thickness is initially 1 surface atomic layer.
#   - During production, the center heat source is re-applied every
#     HeatingChunkSteps. This mimics continuous laser heating rather than a
#     single initial pulse.
#
# Notes:
#   Forcite is classical MD. This script does not model optical absorption,
#   electrons, plasma, or two-temperature physics. It implements a local
#   thermal source in the atomic velocities.
# -------------------------------------------------------------------------

# Empty means use the active document from the User menu.
my $InputDocumentName = "C1.xsd";

my $MakeWorkingCopy = "Yes";
my $OutputDocumentName = "";
my $MakeP1BeforeRun = "No";
my $MaxAutoNameLength = 15;      # Includes extension, e.g. W123456.xsd = 11
my $RunId = int(rand(90)) + 10;  # Two digits keep names short but distinct

# Geometry of the surface.
my $SurfaceAxis = "Z";           # Usually Z for a slab with vacuum above it
my $LaserFrom = "Top";           # Top or Bottom side along SurfaceAxis
my $FixedFrom = "Bottom";        # Fixed support side
my $FixedLayerCount = 5;         # Fixed bottom layers
my $HotLayerCount = 3;           # Initial hot-region thickness in atomic layers
my $HotRadius = 6;             # Angstrom
my $LayerTolerance = 0.45;       # Angstrom; increase if layers are rough
my $HeatOnlyElement = "Cu";      # Empty string heats all elements

# If these are undef, the script uses the center of the selected surface layer.
my $CenterA = undef;             # In-plane coordinate 1, e.g. X when SurfaceAxis=Z
my $CenterB = undef;             # In-plane coordinate 2, e.g. Y when SurfaceAxis=Z

# Thermal setup.
my $BaseTemperature = 300.0;     # K for nonfixed atoms outside the hot source
my $HotTemperature = 8000.0;     # K for center heating region
my $RemoveDriftMode = "Lateral"; # None, Lateral, or All
my $RandomSeed = 24681357;

# Forcite settings. For Cu, try SuttonChen or ZhouJohnsonWadley if installed.
my $Forcefield = "SuttonChen";
my $Quality = "Medium";
my $ChargeAssignment = "Forcefield assigned";
my $AssignForcefieldTypes = "Yes";
my $WriteLevel = "Summary";

# Optional structure preparation.
my $RunGeometryOptimization = "No";
my $GeometryOptimizationSteps = 50;

# Optional 300 K pre-equilibration before imposing the local hot source.
my $RunEquilibration = "Yes";
my $EquilibrationSteps = 2000;
my $EquilibrationEnsemble = "NVT";

# Heat-transfer production run.
my $ProductionSteps = 20000;        # Total production MD steps
my $ContinuousHeating = "Yes";      # Yes = reheat center every chunk
my $HeatingChunkSteps = 1000;        # Steps between heat-source re-application
my $CreateCombinedTrajectory = "Yes";
my $CombinedTrajectoryName = "";
my $ProductionEnsemble = "NVE";
my $TimeStep = 1;              # fs. Use 0.25-0.5 fs for strong heating.
my $TrajectoryFrequency = 500;
my $DataFrequency = 100;
my $EnergyDeviation = 1.0e8;     # kcal/mol, high to avoid premature abort

Documents->SaveAllAtEnd = "Yes";
srand($RandomSeed);
my $HasEquilibrated = 0;

my $source_doc = get_input_document($InputDocumentName);
my $doc = prepare_working_document($source_doc, $OutputDocumentName);

if ($MakeP1BeforeRun =~ /^y/i) {
    print "Converting working structure to P1...\n";
    $doc->MakeP1;
}

my $setup = classify_structure($doc);
lock_laser_center($setup);
print_setup_summary($doc, $setup);

apply_fixed_constraints($setup->{fixed_atoms});

configure_forcite_common();

if ($RunGeometryOptimization =~ /^y/i) {
    print "\nRunning constrained geometry optimization...\n";
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

    $setup = classify_structure($doc);
    apply_fixed_constraints($setup->{fixed_atoms});
    $HasEquilibrated = 1;
}

if ($RunEquilibration =~ /^y/i && $EquilibrationSteps > 0) {
    print "\nRunning 300 K pre-equilibration...\n";
    assign_uniform_temperature($setup->{mobile_atoms}, $BaseTemperature);
    Modules->Forcite->Dynamics->Run(
        $doc,
        dynamics_settings($EquilibrationEnsemble,
                          $BaseTemperature,
                          $EquilibrationSteps,
                          $TrajectoryFrequency,
                          $DataFrequency,
                          "Current")
    );

    # Forcite creates a trajectory using the input document name. Save the
    # equilibrated final structure as a new seed to avoid trajectory-name
    # collisions in the production Dynamics run.
    my $production_doc_name = short_document_name("P$RunId", "xsd");
    print "\nSaving equilibrated production seed: $production_doc_name\n";
    $doc = $doc->SaveAs($production_doc_name);

    $setup = classify_structure($doc);
    apply_fixed_constraints($setup->{fixed_atoms});
}

run_heated_production($doc, $setup);
print "Finished center-heating heat-transfer MD.\n";

sub get_input_document {
    my ($name) = @_;

    if (defined($name) && $name ne "") {
        my $doc = eval { $Documents{$name}; };
        die "Could not find project document '$name'. Check the exact name.\n"
            unless defined $doc;
        return $doc;
    }

    my $doc = eval { Documents->ActiveDocument; };
    die "No InputDocumentName was set and no active document was available.\n"
        unless defined $doc;
    return $doc;
}

sub prepare_working_document {
    my ($source, $output_name) = @_;
    return $source unless $MakeWorkingCopy =~ /^y/i;

    if (!defined($output_name) || $output_name eq "") {
        $output_name = short_document_name("W$RunId", "xsd");
    }
    $output_name = compact_document_name($output_name, "xsd");

    print "Creating working copy: $output_name\n";
    return $source->SaveAs($output_name);
}

sub classify_structure {
    my ($doc, $create_sets) = @_;
    $create_sets = "Yes" unless defined $create_sets;
    my @atoms = @{$doc->UnitCell->Atoms};
    die "No atoms were found in the input structure.\n" unless @atoms;

    my @layers = build_layers(\@atoms, $SurfaceAxis, $LayerTolerance);
    die "Only " . scalar(@layers) . " atomic layers were found; cannot fix $FixedLayerCount layers.\n"
        if @layers <= $FixedLayerCount;
    my @fixed_layers = select_side_layers(\@layers, $FixedFrom, $FixedLayerCount);
    my @center_layers = select_side_layers(\@layers, $LaserFrom, 1);

    my @surface_atoms = map { @{$_->{atoms}} } @center_layers;
    my ($plane_a, $plane_b) = plane_axes($SurfaceAxis);
    my ($center_a, $center_b) =
        defined($CenterA) && defined($CenterB)
        ? ($CenterA, $CenterB)
        : plane_center(\@surface_atoms, $plane_a, $plane_b);

    my %fixed_id = map { atom_key($_) => 1 } map { @{$_->{atoms}} } @fixed_layers;

    my (@fixed_atoms, @mobile_atoms);
    foreach my $atom (@atoms) {
        if ($fixed_id{atom_key($atom)}) {
            push @fixed_atoms, $atom;
            next;
        }

        push @mobile_atoms, $atom;
    }

    my ($effective_hot_layer_count, $hot_layers_ref, $hot_atoms_ref) =
        find_hot_atoms_with_adaptive_thickness(
            \@layers,
            \@mobile_atoms,
            $plane_a,
            $plane_b,
            $center_a,
            $center_b
        );

    my @hot_layers = @$hot_layers_ref;
    my @hot_atoms = @$hot_atoms_ref;
    my %hot_id = map { atom_key($_) => 1 } @hot_atoms;
    my @cold_mobile_atoms;
    foreach my $atom (@mobile_atoms) {
        push @cold_mobile_atoms, $atom unless $hot_id{atom_key($atom)};
    }

    die "The fixed layer selection is empty. Check FixedFrom and FixedLayerCount.\n"
        unless @fixed_atoms;
    die "The mobile selection is empty. Reduce FixedLayerCount.\n"
        unless @mobile_atoms;

    if ($create_sets =~ /^y/i) {
        my $label = clean_name($doc->Name) . "_" . time() . "_" . int(rand(1000000));
        create_atom_set($doc, "${label}_FixedBottom${FixedLayerCount}Layers", \@fixed_atoms);
        create_atom_set($doc, "${label}_Mobile300K", \@mobile_atoms);
        create_atom_set($doc, "${label}_HotCenter_R${HotRadius}_L${effective_hot_layer_count}", \@hot_atoms);
    }

    return {
        atoms => \@atoms,
        layers => \@layers,
        fixed_atoms => \@fixed_atoms,
        mobile_atoms => \@mobile_atoms,
        cold_mobile_atoms => \@cold_mobile_atoms,
        hot_atoms => \@hot_atoms,
        fixed_layers => \@fixed_layers,
        hot_layers => \@hot_layers,
        effective_hot_layer_count => $effective_hot_layer_count,
        plane_a => $plane_a,
        plane_b => $plane_b,
        center_a => $center_a,
        center_b => $center_b
    };
}

sub build_layers {
    my ($atoms, $axis, $tol) = @_;
    my @sorted = sort { atom_coord($a, $axis) <=> atom_coord($b, $axis) } @$atoms;
    my @layers;

    foreach my $atom (@sorted) {
        my $coord = atom_coord($atom, $axis);

        if (!@layers || abs($coord - $layers[-1]->{center}) > $tol) {
            push @layers, { center => $coord, atoms => [$atom] };
            next;
        }

        push @{$layers[-1]->{atoms}}, $atom;
        my $n = scalar @{$layers[-1]->{atoms}};
        $layers[-1]->{center} = (($layers[-1]->{center} * ($n - 1)) + $coord) / $n;
    }

    return @layers;
}

sub find_hot_atoms_with_adaptive_thickness {
    my ($layers, $mobile_atoms, $plane_a, $plane_b, $center_a, $center_b) = @_;

    my $start_count = $HotLayerCount;
    $start_count = 1 if $start_count < 1;
    my $max_count = scalar @$layers;

    for (my $count = $start_count; $count <= $max_count; ++$count) {
        my @hot_layers = select_side_layers($layers, $LaserFrom, $count);
        my @layer_atoms = map { @{$_->{atoms}} } @hot_layers;
        my %hot_layer_id = map { atom_key($_) => 1 } @layer_atoms;

        my @hot_atoms;
        foreach my $atom (@$mobile_atoms) {
            next unless $hot_layer_id{atom_key($atom)};
            next unless should_heat_atom($atom, $HeatOnlyElement);
            next unless plane_distance($atom, $plane_a, $plane_b,
                                       $center_a, $center_b) <= $HotRadius;
            push @hot_atoms, $atom;
        }

        if (@hot_atoms) {
            if ($count > $HotLayerCount) {
                print "Hot center was empty at $HotLayerCount layer(s); ";
                print "using $count layer(s) instead.\n";
            }
            return ($count, \@hot_layers, \@hot_atoms);
        }
    }

    die "The hot center selection is empty even after expanding from "
      . "$HotLayerCount layer(s) to $max_count layer(s). Increase HotRadius, "
      . "check SurfaceAxis/LaserFrom, or check whether the chosen element exists.\n";
}

sub select_side_layers {
    my ($layers, $side, $count) = @_;
    my @selected;

    if ($side =~ /^top$/i) {
        for (my $i = @$layers - $count; $i < @$layers; ++$i) {
            push @selected, $layers->[$i];
        }
        return @selected;
    }

    if ($side =~ /^bottom$/i) {
        for (my $i = 0; $i < $count; ++$i) {
            push @selected, $layers->[$i];
        }
        return @selected;
    }

    die "Layer side must be Top or Bottom.\n";
}

sub print_setup_summary {
    my ($doc, $setup) = @_;
    my @layers = @{$setup->{layers}};

    print "\nCenter-heating setup\n";
    print "Working document: " . $doc->Name . "\n";
    print "Total atoms: " . scalar(@{$setup->{atoms}}) . "\n";
    print "Detected atomic layers along $SurfaceAxis: " . scalar(@layers) . "\n";
    print "Layer tolerance: $LayerTolerance Angstrom\n";
    print "Fixed atoms: " . scalar(@{$setup->{fixed_atoms}}) . "\n";
    print "Mobile 300 K atoms: " . scalar(@{$setup->{mobile_atoms}}) . "\n";
    print "Hot center atoms: " . scalar(@{$setup->{hot_atoms}}) . "\n";
    print "Surface center in plane " . $setup->{plane_a} . $setup->{plane_b};
    print ": (" . $setup->{center_a} . ", " . $setup->{center_b} . ")\n";
    print "Hot radius: $HotRadius Angstrom\n";
    print "Requested hot layer count: $HotLayerCount\n";
    print "Effective hot layer count: " . $setup->{effective_hot_layer_count} . "\n";
    print "Fixed layer count: $FixedLayerCount\n";
    print "Base temperature: $BaseTemperature K\n";
    print "Hot temperature: $HotTemperature K\n";
    print "Forcefield: $Forcefield\n";
}

sub lock_laser_center {
    my ($setup) = @_;
    if (!defined($CenterA) || !defined($CenterB)) {
        $CenterA = $setup->{center_a};
        $CenterB = $setup->{center_b};
        print "Locked laser spot center at ($CenterA, $CenterB) in the surface plane.\n";
    }
}

sub apply_fixed_constraints {
    my ($fixed_atoms) = @_;
    foreach my $atom (@$fixed_atoms) {
        $atom->Velocity = Point(X => 0.0, Y => 0.0, Z => 0.0);
        $atom->Fix("XYZ");
    }
}

sub assign_initial_temperature {
    my ($mobile_atoms, $hot_atoms, $base_temp, $hot_temp) = @_;

    my %hot_id = map { atom_key($_) => 1 } @$hot_atoms;

    foreach my $atom (@$mobile_atoms) {
        my $temperature = $hot_id{atom_key($atom)} ? $hot_temp : $base_temp;
        my ($vx, $vy, $vz) = random_velocity_vector($atom, $temperature);
        $atom->Velocity = Point(X => $vx, Y => $vy, Z => $vz);
    }

    remove_center_of_mass_drift($mobile_atoms, $RemoveDriftMode, $SurfaceAxis);
}

sub assign_uniform_temperature {
    my ($atoms, $temperature) = @_;

    foreach my $atom (@$atoms) {
        my ($vx, $vy, $vz) = random_velocity_vector($atom, $temperature);
        $atom->Velocity = Point(X => $vx, Y => $vy, Z => $vz);
    }

    remove_center_of_mass_drift($atoms, $RemoveDriftMode, $SurfaceAxis);
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

sub dynamics_settings {
    my ($ensemble, $temperature, $steps, $traj_freq, $data_freq, $initial_velocities) = @_;

    return Settings(
        CurrentForcefield => $Forcefield,
        Quality => $Quality,
        WriteLevel => $WriteLevel,
        AssignForcefieldTypes => $AssignForcefieldTypes,
        ChargeAssignment => $ChargeAssignment,
        Ensemble3D => $ensemble,
        Ensemble2D => $ensemble,
        Ensemble0D => $ensemble,
        Temperature => $temperature,
        InitialVelocities => $initial_velocities,
        RandomSeed => $RandomSeed,
        TimeStep => $TimeStep,
        NumberOfSteps => $steps,
        TrajectoryFrequency => $traj_freq,
        DataFrequency => $data_freq,
        WriteVelocities => "Yes",
        WriteForces => "Yes",
        WriteTemperatureData => "Yes",
        WriteEnergyData => "Yes",
        EnergyDeviation => $EnergyDeviation,
        TrajectoryRestart => "No",
        AppendTrajectory => "No"
    );
}

sub run_heated_production {
    my ($doc, $setup) = @_;

    if (!$HasEquilibrated) {
        print "\nNo pre-equilibration was run; assigning $BaseTemperature K to all mobile atoms.\n";
        assign_uniform_temperature($setup->{mobile_atoms}, $BaseTemperature);
    }

    my $chunk_steps = ($ContinuousHeating =~ /^y/i) ? $HeatingChunkSteps : $ProductionSteps;
    die "HeatingChunkSteps must be >= 1.\n" if $chunk_steps < 1;
    die "ProductionSteps must be >= 1.\n" if $ProductionSteps < 1;

    my $combined_doc = undef;
    if ($CreateCombinedTrajectory =~ /^y/i) {
        my $combined_name = $CombinedTrajectoryName;
        $combined_name = combined_trajectory_name($ProductionSteps, $chunk_steps)
            if !defined($combined_name) || $combined_name eq "";
        $combined_name = compact_document_name($combined_name, "xtd");
        print "\nCreating combined production trajectory: $combined_name\n";
        $combined_doc = Documents->New($combined_name);
    }

    print "\nRunning continuous center-heating production dynamics...\n";
    print "Total steps: $ProductionSteps; heating interval: $chunk_steps steps\n";

    my $remaining_steps = $ProductionSteps;
    my $cycle = 1;
    my $last_trajectory = undef;

    while ($remaining_steps > 0) {
        my $steps_this_cycle = $remaining_steps < $chunk_steps
            ? $remaining_steps
            : $chunk_steps;
        my $traj_freq = safe_frequency($steps_this_cycle, $TrajectoryFrequency);
        my $data_freq = safe_frequency($steps_this_cycle, $DataFrequency);

        $setup = classify_structure($doc, "No");
        apply_fixed_constraints($setup->{fixed_atoms});

        print "\nHeating cycle $cycle: applying hot source to ";
        print scalar(@{$setup->{hot_atoms}}) . " atoms ";
        print "(" . $setup->{effective_hot_layer_count} . " layer(s)) ";
        print "at $HotTemperature K; ";
        print "running $steps_this_cycle NVE steps.\n";
        apply_hot_source($setup->{hot_atoms}, $HotTemperature);
        remove_center_of_mass_drift($setup->{mobile_atoms},
                                    $RemoveDriftMode,
                                    $SurfaceAxis);

        my $results = Modules->Forcite->Dynamics->Run(
            $doc,
            dynamics_settings($ProductionEnsemble,
                              $BaseTemperature,
                              $steps_this_cycle,
                              $traj_freq,
                              $data_freq,
                              "Current")
        );

        $last_trajectory = eval { $results->Trajectory; };
        if (defined $last_trajectory) {
            print "Cycle trajectory: " . $last_trajectory->Name . "\n";
            append_to_combined_trajectory($combined_doc, $last_trajectory)
                if defined $combined_doc;
        }

        $remaining_steps -= $steps_this_cycle;
        last if $remaining_steps <= 0;

        my $next_seed_name = cycle_seed_name($cycle + 1);
        print "Saving next cycle seed: $next_seed_name\n";
        $doc = $doc->SaveAs($next_seed_name);
        ++$cycle;
    }

    if (defined $combined_doc) {
        $combined_doc->Save;
        print "Combined trajectory saved as: " . $combined_doc->Name . "\n";
    }
    print "Last production trajectory: " . $last_trajectory->Name . "\n"
        if defined $last_trajectory;
}

sub apply_hot_source {
    my ($hot_atoms, $temperature) = @_;
    die "Hot atom selection is empty; cannot apply continuous heating.\n"
        unless @$hot_atoms;

    foreach my $atom (@$hot_atoms) {
        my ($vx, $vy, $vz) = random_velocity_vector($atom, $temperature);
        $atom->Velocity = Point(X => $vx, Y => $vy, Z => $vz);
    }
}

sub safe_frequency {
    my ($steps, $requested) = @_;
    return $steps if !defined($requested) || $requested < 1 || $requested > $steps;
    return $requested;
}

sub append_to_combined_trajectory {
    my ($combined_doc, $chunk_trajectory) = @_;
    return unless defined $combined_doc && defined $chunk_trajectory;

    my $ok = eval {
        $combined_doc->Trajectory->AppendFramesFrom($chunk_trajectory);
        1;
    };
    print "Warning: could not append chunk trajectory to combined trajectory: $@\n"
        unless $ok;
}

sub random_velocity_vector {
    my ($atom, $temperature) = @_;

    my $mass = $atom->Mass;
    die "Atom mass is zero or undefined for " . $atom->ElementSymbol . ".\n"
        unless defined($mass) && $mass > 0.0;

    # sigma^2 = kB*T/(m*C), C converts amu*(Angstrom/ps)^2 to kcal/mol.
    my $kb = 0.00198720425864083;
    my $conv = 0.00239005736153349;
    my $sigma = sqrt($kb * $temperature / ($conv * $mass));

    return (gaussian() * $sigma,
            gaussian() * $sigma,
            gaussian() * $sigma);
}

sub gaussian {
    my $u1 = rand();
    my $u2 = rand();
    $u1 = 1.0e-12 if $u1 < 1.0e-12;
    return sqrt(-2.0 * log($u1)) * cos(2.0 * 3.141592653589793 * $u2);
}

sub remove_center_of_mass_drift {
    my ($atoms, $mode, $surface_axis) = @_;
    return if !defined($mode) || $mode =~ /^none$/i;
    return unless @$atoms;

    my ($m_total, $px, $py, $pz) = (0.0, 0.0, 0.0, 0.0);
    foreach my $atom (@$atoms) {
        my $v = $atom->Velocity;
        next unless defined $v;
        my $m = $atom->Mass;
        $m_total += $m;
        $px += $m * $v->X;
        $py += $m * $v->Y;
        $pz += $m * $v->Z;
    }
    return if $m_total <= 0.0;

    my ($dvx, $dvy, $dvz) = ($px / $m_total, $py / $m_total, $pz / $m_total);
    foreach my $atom (@$atoms) {
        my $v = $atom->Velocity;
        next unless defined $v;
        my ($vx, $vy, $vz) = ($v->X, $v->Y, $v->Z);

        if ($mode =~ /^all$/i) {
            $vx -= $dvx;
            $vy -= $dvy;
            $vz -= $dvz;
        } elsif ($mode =~ /^lateral$/i) {
            if ($surface_axis =~ /^x$/i) {
                $vy -= $dvy;
                $vz -= $dvz;
            } elsif ($surface_axis =~ /^y$/i) {
                $vx -= $dvx;
                $vz -= $dvz;
            } else {
                $vx -= $dvx;
                $vy -= $dvy;
            }
        } else {
            die "RemoveDriftMode must be None, Lateral, or All.\n";
        }

        $atom->Velocity = Point(X => $vx, Y => $vy, Z => $vz);
    }
}

sub atom_coord {
    my ($atom, $axis) = @_;
    return $atom->X if $axis =~ /^x$/i;
    return $atom->Y if $axis =~ /^y$/i;
    return $atom->Z if $axis =~ /^z$/i;
    die "SurfaceAxis must be X, Y, or Z.\n";
}

sub atom_plane_coord {
    my ($atom, $axis) = @_;
    return $atom->X if $axis =~ /^x$/i;
    return $atom->Y if $axis =~ /^y$/i;
    return $atom->Z if $axis =~ /^z$/i;
    die "Plane axis must be X, Y, or Z.\n";
}

sub plane_axes {
    my ($surface_axis) = @_;
    return ("Y", "Z") if $surface_axis =~ /^x$/i;
    return ("X", "Z") if $surface_axis =~ /^y$/i;
    return ("X", "Y") if $surface_axis =~ /^z$/i;
    die "SurfaceAxis must be X, Y, or Z.\n";
}

sub plane_center {
    my ($atoms, $axis_a, $axis_b) = @_;
    die "Cannot determine surface center from an empty atom set.\n" unless @$atoms;

    my $min_a = atom_plane_coord($atoms->[0], $axis_a);
    my $max_a = $min_a;
    my $min_b = atom_plane_coord($atoms->[0], $axis_b);
    my $max_b = $min_b;

    foreach my $atom (@$atoms) {
        my $a = atom_plane_coord($atom, $axis_a);
        my $b = atom_plane_coord($atom, $axis_b);
        $min_a = $a if $a < $min_a;
        $max_a = $a if $a > $max_a;
        $min_b = $b if $b < $min_b;
        $max_b = $b if $b > $max_b;
    }

    return (($min_a + $max_a) / 2.0, ($min_b + $max_b) / 2.0);
}

sub plane_distance {
    my ($atom, $axis_a, $axis_b, $center_a, $center_b) = @_;
    my $da = atom_plane_coord($atom, $axis_a) - $center_a;
    my $db = atom_plane_coord($atom, $axis_b) - $center_b;
    return sqrt($da * $da + $db * $db);
}

sub should_heat_atom {
    my ($atom, $element) = @_;
    return 1 if !defined($element) || $element eq "";
    return $atom->ElementSymbol eq $element;
}

sub create_atom_set {
    my ($doc, $name, $atoms) = @_;
    return unless @$atoms;
    eval { $doc->CreateSet($name, $atoms); };
    print "Warning: could not create set '$name': $@\n" if $@;
}

sub atom_key {
    my ($atom) = @_;
    return "$atom";
}

sub clean_name {
    my ($name) = @_;
    $name =~ s/[^A-Za-z0-9.]+/_/g;
    $name =~ s/^[_\.]+|[_\.]+$//g;
    return $name;
}

sub unique_xsd_name {
    my ($prefix) = @_;
    return short_document_name($prefix, "xsd");
}

sub unique_xtd_name {
    my ($prefix) = @_;
    return short_document_name($prefix, "xtd");
}

sub short_document_name {
    my ($prefix, $extension) = @_;
    $prefix = compact_base_name($prefix || "D");
    $prefix = substr($prefix, 0, $MaxAutoNameLength - length($extension) - 1);
    $prefix = "D" if $prefix eq "";

    my $max_base_len = $MaxAutoNameLength - length($extension) - 1;
    die "MaxAutoNameLength is too small for .$extension files.\n"
        if $max_base_len < 4;

    for (my $attempt = 0; $attempt < 1000; ++$attempt) {
        my $token = $attempt == 0 ? "" : sprintf("%02d", $attempt);
        my $base = substr($prefix . $token, 0, $max_base_len);
        my $name = "$base.$extension";

        next if document_exists($name);
        die "Internal naming error: '$name' exceeds $MaxAutoNameLength characters.\n"
            if length($name) > $MaxAutoNameLength;
        return $name;
    }

    die "Could not create a unique short .$extension document name.\n";
}

sub cycle_seed_name {
    my ($cycle) = @_;
    return short_document_name("C$RunId.$cycle", "xsd");
}

sub combined_trajectory_name {
    my ($total_steps, $chunk_steps) = @_;
    my $cycles = int(($total_steps + $chunk_steps - 1) / $chunk_steps);
    return short_document_name("T1.$cycles", "xtd");
}

sub compact_document_name {
    my ($name, $default_extension) = @_;
    $name = "" unless defined $name;
    $name =~ s/\\/\//g;
    $name =~ s!.*/!!;

    my $extension = $default_extension;
    if ($name =~ s/\.([A-Za-z0-9]+)$//) {
        $extension = lc($1);
    }

    my $base = compact_base_name($name);
    my $max_base_len = $MaxAutoNameLength - length($extension) - 1;
    die "MaxAutoNameLength is too small for .$extension files.\n"
        if $max_base_len < 4;
    $base = substr($base, 0, $max_base_len) if length($base) > $max_base_len;
    $base = "D" if $base eq "";

    my $candidate = "$base.$extension";
    die "Internal naming error: '$candidate' exceeds $MaxAutoNameLength characters.\n"
        if length($candidate) > $MaxAutoNameLength;
    return $candidate;
}

sub compact_base_name {
    my ($name) = @_;
    $name = clean_name($name);
    $name = compact_sequential_digit_runs($name);
    $name =~ s/_+/_/g;
    $name =~ s/^_+|_+$//g;
    return $name;
}

sub compact_sequential_digit_runs {
    my ($text) = @_;
    $text =~ s/(\d{4,})/compact_digit_run($1)/eg;
    return $text;
}

sub compact_digit_run {
    my ($digits) = @_;

    for (my $end = 999; $end >= 2; --$end) {
        my $seq = join("", 1 .. $end);
        return "1.$end" if $seq eq $digits;
    }

    my @chars = split //, $digits;
    if (@chars >= 4) {
        my $ok = 1;
        for (my $i = 0; $i < @chars; ++$i) {
            if ($chars[$i] != (($i + 1) % 10)) {
                $ok = 0;
                last;
            }
        }
        return "1." . scalar(@chars) if $ok;
    }

    return $digits;
}

sub document_exists {
    my ($name) = @_;
    my $doc = eval { $Documents{$name}; };
    return defined $doc ? 1 : 0;
}
