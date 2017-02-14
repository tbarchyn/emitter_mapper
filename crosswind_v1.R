# emitter_mapper
# Copyright 2017 Thomas E. Barchyn
# Contact: Thomas E. Barchyn [tbarchyn@gmail.com]

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# Please familiarize yourself with the license of this tool, available
# in the distribution with the filename: /docs/license.txt
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

library (spatstat)
library (lattice)

#############################################################
# parameters
# model space parameters
dim_dw <- 200
dim_cw <- 100

# flux measurement
survey_dist <- 20
survey_locs <- seq (1, dim_dw, survey_dist)
survey_locs <- c(survey_locs, dim_dw)

#############################################################
# create a fake emissions surface with some random hot spots
vals <- rnorm (dim_dw * dim_cw, mean = 0.2, sd = 1)
vals[vals < 0] <- 0
em <- matrix (vals, nrow = dim_dw, ncol = dim_cw)
em[20:31, 34:45] <- 1.8
em[80:83, 60:63] <- 3.0
em[20:25, 1:3] <- 5.0
em[100:150, 80:83] <- 1.9
em <- blur (as.im(em), sigma = 5, bleed = TRUE)
em <- as.matrix (em)
levelplot (em, 
           panel = function (...) {
               panel.levelplot (...)
               panel.abline (v = survey_locs, col = 'grey')
           },
           col.regions = colorRampPalette(c('blue', 'white', 'red'))(100),
           main = 'emissions flux'
           )

#############################################################
# plot accumulation downwind
emc <- apply (em, cumsum, MARGIN = 2)
levelplot (emc, 
           panel = function (...) {
               panel.levelplot (...)
               panel.abline (v = survey_locs, col = 'grey')
           },
           col.regions = colorRampPalette(c('purple', 'white', 'green'))(100),
           main = 'column averaged concentrations'
           )

#############################################################
# extract flux planes and create measured flux planes
num_planes <- length (survey_locs)
real_flux_planes <- matrix(NA, nrow = num_planes, ncol = dim_cw)
for (i in 1:num_planes) {
    real_flux_planes[i, ] <- emc[survey_locs[i], ]
}

# slightly obsfucate the real flux planes to represent measurement error
meas_flux_planes <- real_flux_planes + rnorm (num_planes * length (real_flux_planes[1,]),
                                              mean = 0, sd = 0.5)
meas_flux_planes[meas_flux_planes < 0] <- 0.0

# plot a synthetic crossflux survey
plot (meas_flux_planes[2, ], col = 'blue', xlab = 'crosswind distance', ylab = 'concentrations')
lines (real_flux_planes[2, ], col = 'red')

#############################################################
# reconstruct area emissions
em_meas <- em * NA
for (j in num_planes:2) {
    enhancement <- meas_flux_planes[j, ] - meas_flux_planes[j-1, ]
    for (i in 1:dim_cw) {
        survey_dist <- survey_locs[j] - survey_locs[j-1]
        em_meas[survey_locs[j]:survey_locs[j-1], i] <- enhancement[i] / survey_dist
    }
}

levelplot (em_meas, 
           panel = function (...) {
               panel.levelplot (...)
               panel.abline (v = survey_locs, col = 'grey')
           },
           col.regions = colorRampPalette(c('blue', 'white', 'red'))(100),
           main = 'reconstructed emissions flux'
)













