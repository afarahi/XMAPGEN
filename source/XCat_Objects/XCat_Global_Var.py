from XCat_Utilities import read_data_string, read_data_int 

FDIR_HALOS = read_data_string(tag_name='Halo_dir', file_name='./Directories.xml')
REALIZATION_ID = read_data_int(tag_name='Realization_id', file_name='./Directories.xml')
FDIR_SB_MAP_SAVE = read_data_string(tag_name='SB_Map_dir', file_name='./Directories.xml')
FDIR_SB_MAP_SAVE = FDIR_SB_MAP_SAVE + '%i/'%REALIZATION_ID

FDIR_EVENT_MAP = read_data_string(tag_name='Event_Map_dir', file_name='./Directories.xml')
FDIR_EVENT_MAP = FDIR_EVENT_MAP + '%i/'%REALIZATION_ID



