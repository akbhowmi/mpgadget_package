import numpy
import bigfile
import os

def get_snapshot_redshift_correspondence(output_path):
    output_file_names=os.listdir(output_path)
#    print(output_file_names)
#    print(output_path)
    snapshot_space=[]
    redshift_space=[]
    for name in output_file_names:
        if ('PIG' in name):
          try:
            snapshot_number=name[4:]
            snapshot_space.append(snapshot_number)
          except:
            print("Warning: Ignoring filename:%s"%name)
    snapshot_space=numpy.sort(numpy.array(snapshot_space))
    for snapshot_number in snapshot_space:
            pig=bigfile.BigFile(output_path+'PIG_%s'%snapshot_number)
            header=pig.open('/Header')
            redshift=1./header.attrs['Time'][0]-1 
            redshift_space.append(redshift)   
    return numpy.array(snapshot_space),numpy.array(redshift_space) 


def desired_redshift_to_output_redshift(output_path,desired_redshift,list_all=True,file_format='fof_subfind'): 
    snapshot_space,redshift_space=get_snapshot_redshift_correspondence(output_path)
    redshift_difference=numpy.abs(redshift_space-desired_redshift)
    min_redshift_difference=numpy.amin(redshift_difference)
    output_snapshot=snapshot_space[redshift_difference==min_redshift_difference][0]
    output_redshift=redshift_space[redshift_difference==min_redshift_difference][0]
    if (list_all):
        print("Desired redshift: ",desired_redshift)
        print("Output redshift: ",output_redshift)
        print("Output snapshot: ",output_snapshot)            
    return output_redshift,output_snapshot 



def load_snapshot_header(output_path,desired_redshift):
        output_redshift,snapshot_number=desired_redshift_to_output_redshift(output_path,desired_redshift)
        pig=bigfile.BigFile(output_path+'PIG_%s'%snapshot_number)
        header=pig.open('/Header')
        return header
       
    
def get_box_size(output_path):
    output_file_names=os.listdir(output_path)
    snapshot_space=[]
    redshift_space=[]
    for name in output_file_names:
        if ('PIG' in name):
            snapshot_number=name[4:]
            pig=bigfile.BigFile(output_path+'PIG_%s'%snapshot_number)
            header=pig.open('/Header')
            return header.attrs['BoxSize'][0]


        
def get_particle_property(output_path,particle_property,p_type,desired_redshift):
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift)
    pig=bigfile.BigFile(output_path+'PIG_%s'%output_snapshot)
    return pig.open('%d/%s'%(p_type,particle_property))[:],output_redshift


def get_group_property(output_path,group_property,desired_redshift):
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift)
    pig=bigfile.BigFile(output_path+'PIG_%s'%output_snapshot)
    return pig.open('FOFGroups/%s'%(group_property))[:],output_redshift

    #return il.snapshot.loadSubset(output_path,output_snapshot,p_type,fields=particle_property),output_redshift

def get_particle_property_within_groups(output_path,particle_property,p_type,desired_redshift,group_index):
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift)
    pig=bigfile.BigFile(output_path+'PIG_%s'%output_snapshot)
    GroupID=pig.open('%d/GroupID'%(p_type))[:]
    finalproperty=pig.open('%d/%s'%(p_type,particle_property))[:]   
    return finalproperty[GroupID==GroupID[group_index]],output_redshift


def get_halo_density_profile(output_path,p_type,desired_redshift_of_selected_halo,index_of_selected_halo,min_edge,max_edge,Nbins,CENTER_AROUND='POTENTIAL_MINIMUM',p_id=0):
    from kdcount import correlate
    def min_dis(median_position, position,box_size):
        pos_1=position-median_position
        pos_2=position-median_position+boxsize
        pos_3=position-median_position-boxsize
        new_position_options=numpy.array([pos_1,pos_2,pos_3])
        get_minimum_distance=numpy.argmin(numpy.abs(new_position_options))
        return new_position_options[get_minimum_distance]
    boxsize=get_box_size(output_path)
    particle_property='Position'
    group_positions,output_redshift=get_particle_property_within_groups(output_path,particle_property,p_type,desired_redshift_of_selected_halo,index_of_selected_halo)

    particle_property='Mass'
    group_mass,output_redshift=get_particle_property_within_groups(output_path,particle_property,p_type,desired_redshift_of_selected_halo,index_of_selected_halo) 
    particle_property='Potential'
    group_potential,output_redshift=get_particle_property_within_groups(output_path,particle_property,p_type,desired_redshift_of_selected_halo,index_of_selected_halo)
    if (CENTER_AROUND=='MOST_MASSIVE_BLACKHOLE'):

        particle_property='ID'
        
        bh_IDs,output_redshift=get_particle_property_within_groups(output_path,particle_property,5,desired_redshift_of_selected_halo,index_of_selected_halo)        

        
        
        particle_property='Position'
        
        bh_positions,output_redshift=get_particle_property_within_groups(output_path,particle_property,5,desired_redshift_of_selected_halo,index_of_selected_halo)        
        particle_property='BlackholeMass'

        bh_masses,output_redshift=get_particle_property_within_groups(output_path,particle_property,5,desired_redshift_of_selected_halo,index_of_selected_halo) 
        print("Calculating density around BH with ID:",(bh_IDs[bh_masses==numpy.amax(bh_masses)])[0])
        center=(bh_positions[bh_masses==numpy.amax(bh_masses)])[0]        
    if (CENTER_AROUND=='POTENTIAL_MINIMUM'):
        center=(group_positions[group_potential==numpy.amin(group_potential)])[0]
    transposed_group_positions=numpy.transpose(group_positions)
    vectorized_min_dis = numpy.vectorize(min_dis)
    x_dis=vectorized_min_dis(center[0],transposed_group_positions[0],boxsize)
    y_dis=vectorized_min_dis(center[1],transposed_group_positions[1],boxsize)
    z_dis=vectorized_min_dis(center[2],transposed_group_positions[2],boxsize)
    log_distances=numpy.log10(numpy.sqrt(x_dis**2+y_dis**2+z_dis**2))

    log_distance_bins=numpy.linspace(min_edge,max_edge,Nbins)
    binning=correlate.RBinning(log_distance_bins)
    bin_edges=binning.edges
    bin_centers=binning.centers
    mass_distribution=[]
    for i in range(0,len(bin_edges)-1):
        left=bin_edges[i]
        right=bin_edges[i+1]
        mask=(log_distances>left)&(log_distances<right)
        mass_inside_bin=numpy.sum(group_mass[mask])
        mass_distribution.append(mass_inside_bin)

    mass_distribution=numpy.array(mass_distribution)
    mass_density=mass_distribution/4./3.14/(10**bin_centers)**3/((numpy.diff(bin_centers))[0])/numpy.log(10)
    return bin_centers,mass_distribution,mass_density





