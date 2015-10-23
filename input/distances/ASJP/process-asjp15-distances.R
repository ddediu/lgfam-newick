###############################################################################################################
#
# Convert the distances matrix as created by asjp62 (from listss15.txt) into a proper CSV (TAB-separated)
# matrix distance between languages, using ISO language codes.
#
# Copyright (C) 2013-2014  Dan Dediu
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

library(stringr); # trim spaces from strings


# Output formats:
do.save.as.csv <- FALSE; # save the resulting distances matrix as a CSV (TAB-separated) file (this is really big)
do.save.as.RData <- TRUE; # save the resulting distances matrix as a compressed (xz) RData object


# The log file:
cat( "Converting the ASJP15 distances to a distances matrix.\n", file="./conversion.log", append=FALSE );


# Parse the original file listss15.txt to recover the language info (name + classification, ISO and WALS codes):
lgs <- data.frame( "Language"=rep(NA,100000), "WALS"=NA, "ISO"=NA, stringsAsFactors=FALSE );
con <- file( "./listss15.txt", "r", blocking = FALSE ); # WARNING: first decompress listss15.txt.tar.xz into the same folder as this script! 
k <- 1;
while( length( cur.line <- readLines( con, n=1 ) ) > 0 )
{
  if( length( grep( "{", cur.line, fixed=TRUE, useBytes=TRUE ) ) == 1 )
  {
    # Language line, look at the next line for more info:
    if( length( next.line <- readLines( con, n=1 ) ) == 0 )
    {
      stop( "Error parsing file!\n" );
    }
    # Store this language:
    #lgs <- rbind( lgs,  data.frame( "Language"=cur.line, "WALS"=substr( next.line, 34, 36 ), "ISO"=substr( next.line, 40, 42 ), stringsAsFactors=FALSE ) );
    lgs$Language[k] <- cur.line;
    lgs$WALS[k] <- substr( next.line, 34, 36 );
    lgs$ISO[k]  <- substr( next.line, 40, 42 );
    k <- k+1;
  }
}
close( con );

# Keep only those that have at least WALS or ISO:
lgs <- lgs[ !is.na(lgs$WALS) | !is.na(lgs$ISO), ]; lgs <- lgs[ str_trim(lgs$WALS) != "" | str_trim(lgs$ISO) != "", ];

# Load the matrix (faster processing):
con <- file( "./asjp15-dists.txt", "r", blocking = FALSE );
tmp <- readLines( con, 5 ); # ignore the first 5 lines
asjp15.d <- readLines( con );
close( con );
asjp15.d <- asjp15.d[ -c( length(asjp15.d), length(asjp15.d)-1 ) ]; # get rid of the last two lines:

# Build the corresponding numeric matrix:
asjp15.dm <- matrix( 0.0, nrow=length(asjp15.d), ncol=length(asjp15.d) );
lgs.to.keep <- rep( NA, length(asjp15.d) );

# Each line is structured as SPACE + 26 chars lg info + numbers
# Match the language names and use them as the corresponding row and colnames:
for( i in 1:length(asjp15.d) )
{
  if( i %% 100 == 0 )
  {
    cat( "Entry ", i, " of ", length(asjp15.d), " (", 100*i/length(asjp15.d), "%)...\n" );
  }
  
  lginfo <- substr( asjp15.d[i], 2, 27 );
  lgpos  <- grep( lginfo, lgs$Language, fixed=TRUE, useBytes=TRUE );
  if( length(lgpos) < 1 )
  {
    cat( "\tLanguage ", lginfo, " was not found; skipping...\n", file="./conversion.log", append=TRUE );
  } else
  {
    if( length(lgpos) > 1 )
    {
      # Pick the first one:
      cat( "\tFound duplicated info for ", lginfo, ": pick the first entry...\n", file="./conversion.log", append=TRUE );
      lgpos <- lgpos[1];
    } 
    
    if( lgs$ISO[lgpos] != "" )
    {
      # Keep this entry (and store the index of the language):
      lgs.to.keep[i] <- lgpos;
      
      # Parse the actual distances and store them:
      dists <- str_trim(substr( asjp15.d[i], 28, nchar(asjp15.d[i]) ));
      dists <- as.numeric( strsplit( dists, "[[:blank:]]+" )[[1]] );
      asjp15.dm[i,1:length(dists)] <- dists;
    }
  }
}
# Remove languages without ISO:
asjp15.dm <- asjp15.dm[ !is.na(lgs.to.keep), !is.na(lgs.to.keep) ]; lgs.to.keep <- lgs.to.keep[ !is.na(lgs.to.keep) ];
# Remove duplicated languages:
asjp15.dm <- asjp15.dm[ !duplicated(lgs.to.keep), !duplicated(lgs.to.keep) ]; lgs.to.keep <- lgs.to.keep[ !duplicated(lgs.to.keep) ];
# Remove languages with the same ISO code:
asjp15.dm <- asjp15.dm[ !duplicated(lgs$ISO[lgs.to.keep]), !duplicated(lgs$ISO[lgs.to.keep]) ]; lgs.to.keep <- lgs.to.keep[ !duplicated(lgs$ISO[lgs.to.keep]) ];
# Set the row and col names to the language ISO:
rownames(asjp15.dm) <- (colnames(asjp15.dm) <- lgs$ISO[lgs.to.keep]);
# Make the matrix a full distances matrix:
asjp15.dm <- asjp15.dm + t(asjp15.dm);
# And normalize it between 0 and 1:
asjp15.dm.min <- min(as.numeric(asjp15.dm)); asjp15.dm.max <- max(as.numeric(asjp15.dm)); asjp15.dm <- (asjp15.dm - asjp15.dm.min)/(asjp15.dm.max - asjp15.dm.min);

if( do.save.as.csv )
{
  # Save this matrix as a CSV file:
  cat( "Writing as CSV (TAB-separated) file)...\n", file="./conversion.log", append=TRUE );
  write.table( asjp15.dm, "./asjp15-dists.csv", quote=FALSE, sep="\t" );
}

if( do.save.as.RData )
{
  # Save as an R object:
  cat( "Writing as R object...\n", file="./conversion.log", append=TRUE );
  save( asjp15.dm, file = "./asjp15-dists.RData",compress="xz" );
}

cat( "Conversion DONE!\n\n", file="./conversion.log", append=TRUE );



