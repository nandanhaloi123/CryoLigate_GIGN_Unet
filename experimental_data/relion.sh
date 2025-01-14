relion_image_handler --i 8xv2_Ligand_Map_filtered_gaus_0.7.mrc --new_box 48 --o 8xv2_Ligand_Map_filtered_gaus_0.7_newbox.mrc

relion_image_handler --i 8xv2_Ligand_Map_filtered_gaus_0.7_newbox.mrc --o 8xv2_Ligand_Map_filtered_gaus_0.7_newscale.mrc --rescale_angpix 0.5

relion_image_handler --i 8xv2_Ligand_Map_filtered_gaus_0.7_newscale.mrc --new_box 48 --o 8xv2_Ligand_Map_filtered_gaus_0.7_newscale_newbox.mrc

