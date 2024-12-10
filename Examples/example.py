from ovito.io import import_file, export_file 
from CalculateHydrogenBonds import CalculateHydrogenBonds

pipeline = import_file("input_file.lammpstrj")  
total_frame = pipeline.source.num_frames
pipeline.modifiers.append(CalculateHydrogenBonds(donor= 5, accep=3,  H_atom=6, da_rcut=3.5, dha_rangle=150))
pipeline.modifiers.clear
export_file(pipeline, "output_file.dat", format = "txt/attr", columns = ["SourceFrame", "Hbond_count", "donor_count"],
                multiple_frames = True, start_frame = 0, end_frame = total_frame, every_nth_frame = 1)