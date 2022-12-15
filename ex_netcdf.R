# library(coldpool)
library(ncdf4)
library(RNetCDF)
# R Studio tip: Text wrapping

# Load the raster brick by accessing it from the coldpool package or loading it from the rds file
ebs_bt_brick <- readRDS(file = "bottom_temperature.rds")
ebs_bt_brick <- coldpool::nbs_ebs_bottom_temperature

plot(ebs_bt_brick, 
     zlim = c(-2, 15)) # Setting a scale for the z dimension of the individual layers

# Raster and raster bricks/stacks use dimension indexing to efficiently structure multidimensional data
# Raster dimensions are defined by an extent, cell dimensions, and a coordinate reference system.
ebs_bt_brick

object.size(ebs_bt_brick*1)
object.size(as.data.frame(ebs_bt_brick, xy = TRUE))

# Fill in the z dimension (time)
names(ebs_bt_brick)
ebs_bt_brick@z

# When using the raster or ncdf4 packages, writing variables or dimensions to netCDF files requires that variable and dimensions are interpretable as float, double, integer, or binary types because of the need to be interoperable among systems and programming languages.
# This z is layers of the brick, not the values of individual layers... yes, it's a bit confusing
ebs_bt_brick <- raster::setZ(ebs_bt_brick, 
                             z = as.numeric(
                               gsub(x = names(ebs_bt_brick), 
                                    pattern = "[^0-9.-]", 
                                    replacement = "")
                               )
                             )

plot(ebs_bt_brick, 
     zlim = c(-2, 15))

ebs_bt_brick@z

# Write the raster to a netCDF file
writeRaster(x = ebs_bt_brick, 
            filename = "temperature.nc", 
            overwrite = TRUE, 
            format = "CDF", 
            varname = "BT", 
            varunit = "degC", 
            longname = "Bottom temperature", 
            xname = "Longitude", 
            yname = "Latitude", 
            zname = "Year", 
            zunit = "numeric")

# Open a connection to the netCDF file using ncdf4
con_ncdf4 <- ncdf4::nc_open(filename = "temperature.nc", 
                            write = FALSE)

# Retrieve the temperature variable 
bt_array <- ncdf4::ncvar_get(nc = con_ncdf4,
                 varid = "BT")

ncdf4::ncatt_get(nc = con_ncdf4,
                             varid = "BT",
                             attname = "proj4")

ncdf4::ncatt_get(nc = con_ncdf4,
                 varid = "BT",
                 attname = "units")

dim(bt_array) # Dimensions are 335 x 335 x 6, i.e. x dimension (longitude), y dimension (latitude), z dimension (year) 

bt_2010 <- bt_array[, , 1] # Select year 2010

plot(
  raster(
    bt_2010,
    crs = "EPSG:3338"
    ), 
     zlim = c(-2,15)
  )

# Dimension variables are accessed using the same approach
ncdf4::ncvar_get(nc = con_ncdf4,
                 varid = "Longitude")

ncdf4::ncvar_get(nc = con_ncdf4,
                 varid = "Latitude")

ncdf4::ncvar_get(nc = con_ncdf4,
                 varid = "Year")

# Add a new attribute to a variable; setting up a connection in write mode
con_ncdf4 <- ncdf4::nc_open(filename = "temperature.nc",
                            write = TRUE)
ncatt_put(nc = con_ncdf4,
          varid = "BT",
          attname = "method",
          attval = "Interpolated from survey gear temperature observations using Ordinary Kriging with a Stein's Matern Variogram")


# Add a global attribute
ncatt_put(con_ncdf4,
          varid = 0,
          attname = "citation",
          attval = "Rohan, S.K., Barnett L.A.K., and Charriere, N. In review. Evaluating approaches to estimating mean temperatures and cold pool area from AFSC bottom trawl surveys of the eastern Bering Sea. U.S. Dep. Commer., NOAA Tech. Mem. NMFS-AFSC-456, XX p.")



# Add a new variable with the same dimensions as the original layer (surface temperature) ----
ebs_st_brick <- readRDS(file = "surface_temperature.rds")
# ebs_st_brick <- coldpool::nbs_ebs_surface_temperature

saveRDS(ebs_st_brick, file = "surface_temperature.rds")

plot(ebs_st_brick, zlim = c(-2, 18))

x_dim <- con_ncdf4$dim[['Longitude']]
y_dim <- con_ncdf4$dim[['Latitude']]
t_dim <- con_ncdf4$dim[['Year']]

# Setup a variable by defining the name, units, dimensions, and value for missing data
sst_var <- ncvar_def(name = "SST", 
                   units = "degC", 
                   dim = list(x_dim, y_dim, t_dim),
                   missval = -3.39999995214436e+38)

# Add the variable to the netCDF file in 'define' mode
con_ncdf4 <- ncvar_add(con_ncdf4, sst_var)

# Revert to write mode
con_ncdf4 <- ncdf4::nc_open(filename = "temperature.nc",
                            write = TRUE)

# Convert the SST brick to a 3D array and add it to the netCDF file
ncvar_put(nc = con_ncdf4, 
          varid = "SST", 
          vals = as.array(ebs_st_brick))
ncatt_put(nc = con_ncdf4, 
          varid = "SST", 
          attname = "long_name", 
          attval  = "Sea surface temperature")

st_array <- ncvar_get(nc = con_ncdf4,
                      varid = "SST")

plot(
  raster(st_array[,,1],
         crs = "EPSG:3338")
)

# Add a variable with differnt dimensions (average bottom temperature) ----

# Select data from years with temperature data that extends into the NBS
# mean_temperature <- coldpool::cold_pool_index |>
#   dplyr::filter(YEAR %in% c(2010, 2017, 2018, 2019, 2021, 2022))

mean_temperature <- readRDS("mean_temperature.rds")

# Setup a variable by defining the name, units, dimensions, and value for missing data
mean_bt_var <- ncvar_def(name = "Mean BT", 
                     units = "degC", 
                     dim = list(t_dim),
                     missval = -3.39999995214436e+38)

# Add the variable to the netCDF file in 'define' mode
con_ncdf4 <- ncvar_add(con_ncdf4, mean_bt_var)

con_ncdf4 <- ncdf4::nc_open(filename = "temperature.nc",
                            write = TRUE)

ncvar_put(nc = con_ncdf4, 
          varid = "Mean BT", 
          vals = mean_temperature$MEAN_GEAR_TEMPERATURE)

ncatt_put(nc = con_ncdf4, 
          varid = "Mean BT", 
          attname = "long_name", 
          attval  = "Mean EBS shelf bottom temperature")

ncvar_get(nc = con_ncdf4, 
          varid = "Mean BT")

nc_close(con_ncdf4)


# The RNetCDF package is similar to ncdf4 but requires users to set variable types when writing variables, attributes, and dimensions. There are slight differnces in the interfaces differ a bit (example below). An advantage of the RNetCDF package is that it accomodates more variable types, including characters. This can be useful for representing dates and times since ncdf4 can only represent time in relation to a reference point (a geoscience standard practice) since it only supports floats, integers, and binary variables. 


con_rnetcdf <- RNetCDF::open.nc(con = "temperature.nc")
con_rnetcdf
print.nc(con_rnetcdf)

var.get.nc(ncfile = con_rnetcdf, variable = "SST")
att.get.nc(ncfile = con_rnetcdf, variable = "SST", attribute = "units")

# The goal for this demo was to provide an intro to the structure and logic of multidimensional data formats. So, we kinda 'cheated' by using the raster package to setup the netCDF file. Setting up a netCDF file from scratch requires defining dimensions using nc_create() and ncdim_def(). The approach to using these functions is similar as creating variables since dimensions are essentially just variables.