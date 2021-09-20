def add_law_beta_drives(net, beta_start, strength=1.0):
    # Distal Drive
    weights_ampa_d1 = {'L2_basket': 0.00032 * strength, 'L2_pyramidal': 0.00008,
                       'L5_pyramidal': 0.00004}
    syn_delays_d1 = {'L2_basket': 0.5, 'L2_pyramidal': 0.5,
                     'L5_pyramidal': 0.5}
    net.add_bursty_drive(
        'beta_dist', tstart=beta_start, tstart_std=0., tstop=beta_start + 50.,
        burst_rate=1., burst_std=10., numspikes=2, spike_isi=10, n_drive_cells=10,
        location='distal', weights_ampa=weights_ampa_d1,
        synaptic_delays=syn_delays_d1, seedcore=20)

    # Proximal Drive
    weights_ampa_p1 = {'L2_basket': 0.00004 * strength, 'L2_pyramidal': 0.00002,
                       'L5_basket': 0.00002, 'L5_pyramidal': 0.00002}
    syn_delays_p1 = {'L2_basket': 0.1, 'L2_pyramidal': 0.1,
                     'L5_basket': 1.0, 'L5_pyramidal': 1.0}

    net.add_bursty_drive(
        'beta_prox', tstart=beta_start, tstart_std=0., tstop=beta_start + 50.,
        burst_rate=1., burst_std=20., numspikes=2, spike_isi=10, n_drive_cells=10,
        location='proximal', weights_ampa=weights_ampa_p1,
        synaptic_delays=syn_delays_p1, seedcore=20)

    return net


def add_supra_beta_drives(net, beta_start, strength=1.0):
    """Add beta drives using parameters from Jones 2009 ERP drives21"""
    # Add beta evoking drives
    #Distal 
    weights_ampa_d1 = {'L2_basket': 0.006562 * strength, 'L2_pyramidal': 0.000007 * strength,
                    'L5_pyramidal': 0.142300 * strength}
    weights_nmda_d1 = {'L2_basket': 0.019482 * strength, 'L2_pyramidal': 0.004317 * strength,
                    'L5_pyramidal': 0.080074 * strength}
    synaptic_delays_d1 = {'L2_basket': 0.1, 'L2_pyramidal': 0.1,
                        'L5_pyramidal': 0.1}
    net.add_evoked_drive(
        'beta_dist', mu=beta_start, sigma=3.85, numspikes=1, weights_ampa=weights_ampa_d1,
        weights_nmda=weights_nmda_d1, location='distal',
        synaptic_delays=synaptic_delays_d1, seedcore=4)

    # Proximal
    weights_ampa_p2 = {'L2_basket': 0.000003 * strength, 'L2_pyramidal': 1.438840 * strength,
                    'L5_basket': 0.008958 * strength, 'L5_pyramidal': 0.684013 * strength}
    synaptic_delays_prox = {'L2_basket': 0.1, 'L2_pyramidal': 0.1,
                            'L5_basket': 1., 'L5_pyramidal': 1.}
    net.add_evoked_drive(
        'beta_prox', mu=beta_start, sigma=8.33, numspikes=1,
        weights_ampa=weights_ampa_p2, location='proximal',
        synaptic_delays=synaptic_delays_prox, seedcore=4)

def rescale_pyr_morph(net, cell_types, compartment_prop, scaling_factor):
    for cell_type in cell_types:
        for section in net.cell_types[cell_type].sections.values():
            section[compartment_prop] *= scaling_factor