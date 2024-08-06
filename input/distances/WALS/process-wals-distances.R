###############################################################################################################
#
# Compute distances between languages using the WALS database.
#
# Copyright (C) 2014=2024  Dan Dediu
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#################################################################################################################

# Load and merge the wals data in the forma required by the old code:
# iso_code	glottocode	Name	latitude	longitude	genus	family feature_columns_starting_with_X
wals_languages <- read.table( "../../wals/languages.csv", header=TRUE, sep=",", quote="\"", stringsAsFactors=FALSE );
wals_languages <- unique(wals_languages[, c("ID", "ISO639P3code", "Glottocode", "Name", "Latitude", "Longitude", "Genus", "Family")]); # keep only the ones needed
names(wals_languages) <- c("wals_code", "iso_code", "glottocode", "Name", "latitude", "longitude", "genus", "family");
wals_languages <- wals_languages[ !is.na(wals_languages$iso_code) & wals_languages$iso_code != "", ]; # keep only the ones with a ISO code
wals_features <- read.table( "../../wals/parameters.csv", header=TRUE, sep=",", quote="\"", stringsAsFactors=FALSE );
wals_languages <- cbind(wals_languages, matrix(NA, nrow=nrow(wals_languages), ncol=nrow(wals_features))); 
names(wals_languages)[9:ncol(wals_languages)] <- paste0("X", wals_features$ID);
wals_values <- read.table( "../../wals/values.csv", header=TRUE, sep=",", quote="\"", stringsAsFactors=FALSE );
# slow but clear:
for( i in 1:nrow(wals_languages) )
{
  s <- wals_values[ wals_values$Language_ID == wals_languages$wals_code[i], ];
  if( nrow(s) > 0 ) wals_languages[i, paste0("X",s$Parameter_ID) ] <- s$Value;
}


# --- ORIGINAL CODE ---

# Load the wals data:
wals <- wals_languages; # rename it for compatibility with the old code
#wals <- read.table( "../../wals/language.csv", header=TRUE, sep=",", quote="\"", stringsAsFactors=FALSE );

# Keep only the numeric values of the features (discard the explanations):
feat.cols <- grep( "^X", names(wals) );
for( i in 1:nrow(wals) )
{
  #cat( "Processing language ", i, " out of ", nrow(wals), "...\n" );
  for( j in feat.cols )
  {
    if( is.na(wals[i,j]) || wals[i,j] == "" )
    {
      # Mark it as proper NA:
      wals[i,j] <- NA;
    } else
    {
      # Clean the text:
      value <- regexpr( "[[:digit:]]+", as.character(wals[i,j]) );
      if( value[1] == -1 || length(value) == 0 )
      {
        cat( "Malformed value \"", as.character(wals[i,j]), "\"\n" );
        wals[i,j] <- NA;
      } else
      {
        wals[i,j] <- substr( as.character(wals[i,j]), value[1], value[1] + attr(value,"match.length")[1] - 1 );
      }
    }
  }
}
# The wals codes:
wals.codes <- sort( unique( as.character( wals$wals_code ) ) );


# Various types of distances implemented by daisy in cluster:
library(cluster);

# Extract the actual values:
vals <- wals[ , feat.cols ]; rownames(vals) <- wals$wals_code;
vals <- vals[ rownames(vals) %in% wals.codes, ];
vals <- vals[ order( rownames(vals) ), ];

# For gower make sure they are factors:
for( i in 1:ncol(vals) )
{
  vals[,i] <- as.factor(vals[,i]);
}
wals.gower.dm <- daisy( vals, metric="gower", stand=TRUE );

# For manhattan and euclidean they are numeric:
for( i in 1:ncol(vals) )
{
  vals[,i] <- as.numeric(as.character(vals[,i]));
}
wals.manhattan.dm <- daisy( vals, metric="manhattan", stand=FALSE );
wals.euclidean.dm <- daisy( vals, metric="euclidean", stand=TRUE );


# Missing data imputation using the mode:
vals <- wals[ , feat.cols ]; rownames(vals) <- wals$wals_code;
vals <- vals[ rownames(vals) %in% wals.codes, ];
vals <- vals[ order( rownames(vals) ), ];
for( i in 1:ncol(vals) )
{
  # get the mode:
  ux <- unique(vals[,i]); ux <- ux[ !is.na(ux) ];
  if( length(ux) > 0 )
  {
    # Good, some non-missing data:
    ux <- ux[ which.max( tabulate( match(vals[,i], ux) ) ) ];
  }
  vals[is.na(vals[,i]),i] <- ux;
}

# For gower make sure they are factors:
for( i in 1:ncol(vals) )
{
  vals[,i] <- as.factor(vals[,i]);
}
wals.gower.mode.dm <- daisy( vals, metric="gower", stand=TRUE );

# For manhattan and euclidean they are numeric:
for( i in 1:ncol(vals) )
{
  vals[,i] <- as.numeric(as.character(vals[,i]));
}
wals.manhattan.mode.dm <- daisy( vals, metric="manhattan", stand=FALSE );
wals.euclidean.mode.dm <- daisy( vals, metric="euclidean", stand=TRUE );


# Mantel correlations between these distances:
if( FALSE ) # takes some CPU time to run!
{
  library(vegan);
  # With missing data:
  mantel( wals.gower.dm,     wals.manhattan.dm, permutations=999, na.rm=TRUE, parallel=6 ); # r=0.65, p=0.001  <--- gower and manhattan are very similar
  mantel( wals.gower.dm,     wals.euclidean.dm, permutations=999, na.rm=TRUE, parallel=6 ); # r=0.40, p=0.001
  mantel( wals.manhattan.dm, wals.euclidean.dm, permutations=999, na.rm=TRUE, parallel=6 ); # r=0.41, p=0.001
  
  # With mode imputation:
  mantel( wals.gower.mode.dm,     wals.manhattan.mode.dm, permutations=999, na.rm=TRUE, parallel=6 ); # r=0.95, p=0.001  <--- gower and manhattan are very similar
  mantel( wals.gower.mode.dm,     wals.euclidean.mode.dm, permutations=999, na.rm=TRUE, parallel=6 ); # r=0.52, p=0.001
  mantel( wals.manhattan.mode.dm, wals.euclidean.mode.dm, permutations=999, na.rm=TRUE, parallel=6 ); # r=0.54, p=0.001
    
  # Between missing data and mode imputation:
  mantel( wals.gower.dm,     wals.gower.mode.dm,     permutations=999, na.rm=TRUE, parallel=6 ); # r=0.15, p=0.001
  mantel( wals.euclidean.dm, wals.euclidean.mode.dm, permutations=999, na.rm=TRUE, parallel=6 ); # r=0.19, p=0.001
  mantel( wals.manhattan.dm, wals.manhattan.mode.dm, permutations=999, na.rm=TRUE, parallel=6 ); # r=0.27, p=0.001
}


# Save these distances to file:
save( wals.gower.dm,          file="./wals-gower-dm.RData", compress="xz" );
save( wals.manhattan.dm,      file="./wals-manhattan-dm.RData", compress="xz" );
save( wals.euclidean.dm,      file="./wals-euclidean-dm.RData", compress="xz" );
save( wals.gower.mode.dm,     file="./wals-gower-mode-dm.RData", compress="xz" );
save( wals.manhattan.mode.dm, file="./wals-manhattan-mode-dm.RData", compress="xz" );
save( wals.euclidean.mode.dm, file="./wals-euclidean-mode-dm.RData", compress="xz" );





