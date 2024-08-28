// cmass_leaf = cmass_init * (1 - cstorage_max_temp) * ltor / (1 + ltor)
// lai = cmass_leaf * sla
// lai = sla * cmass_init * (1 - cstorage_max_temp) * ltor / (1 + ltor)
const double LAI_MIN = 0.1;
// LAI_MIN = sla * cmass_init * (1 - cstorage_max_temp) * ltor / (1 + ltor)
// cmass_init = LAI_MIN / (sla * (1 - cstorage_max_temp) * ltor / (1 + ltor))
