# mpgadget_package

This package is used to read output files from simulations run using MP-GADGET code (Feng+ 2014).

The package is structured very similar to the pre-existing arepo_package. The function names in both the packages are the same, and they do the exact same job. The only thing that the use needs to modify are the names of different fields depending on how they are stored in MP_GADGET vs AREPO. For example, to the read the positions of ptype = 1, following are the commands in arepo_package vs. mpgadget_package

arepo_package.get_particle_property(output_path,'Coordinates',p_type,desired_redshift)

mpgadget_pacakge.get_particle_property(output_path,'Position',p_type,desired_redshift)

This is because the particle positions are stored as 'Coordinates' in AREPO, and Position' in MP-GADGET


