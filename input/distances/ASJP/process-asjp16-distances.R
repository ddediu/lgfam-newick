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
cat( "Converting the ASJP16 distances to a distances matrix.\n", file="./conversion.log", append=FALSE );


# Parse the original file listss16.txt to recover the language info (name + classification, ISO and WALS codes):
lgs <- data.frame( "Language"=rep(NA,100000), "WALS"=NA, "ISO"=NA, stringsAsFactors=FALSE );
con <- file( "./listss16dd.txt", "r", blocking = FALSE ); # WARNING: first decompress listss16.txt.tar.xz into the same folder as this script! 
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
lgs <- lgs[1:(k-1),];

# Keep only those that have at least WALS or ISO:
for( i in which(str_trim(lgs$WALS) == "" & str_trim(lgs$ISO) == "") )
{
  cat( "\tNo ISO or WALS codes for language ", lgs$Language[i], "; skipping...\n", file="./conversion.log", append=TRUE );
}
lgs <- lgs[ !is.na(lgs$WALS) | !is.na(lgs$ISO), ]; lgs <- lgs[ str_trim(lgs$WALS) != "" | str_trim(lgs$ISO) != "", ];

## filter non-ASCII character in language names:
#library(tools)
#showNonASCII(lgs$Language)

# Load the matrix (faster processing):
con <- file( "./asjp16-dists-dd.txt", "r", blocking = FALSE );
#tmp <- readLines( con, 5 ); # ignore the first 5 lines
asjp16.d <- readLines( con );
close( con );
#asjp16.d <- asjp16.d[ -c( length(asjp16.d), length(asjp16.d)-1 ) ]; # get rid of the last two lines:

# Build the corresponding numeric matrix:
asjp16.dm <- matrix( 0.0, nrow=length(asjp16.d), ncol=length(asjp16.d) );
lgs.to.keep <- rep( NA, length(asjp16.d) );

# Each line is structured as SPACE + 26 chars lg info + numbers
# Match the language names and use them as the corresponding row and colnames:
for( i in 1:length(asjp16.d) )
{
  if( i %% 100 == 0 )
  {
    cat( "Entry ", i, " of ", length(asjp16.d), " (", 100*i/length(asjp16.d), "%)...\n" );
  }
  
  lginfo <- substr( asjp16.d[i], 2, 27 );
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
      dists <- str_trim(substr( asjp16.d[i], 28, nchar(asjp16.d[i]) ));
      dists <- as.numeric( strsplit( dists, "[[:blank:]]+" )[[1]] );
      asjp16.dm[i,1:length(dists)] <- dists;
    }
  }
}
# Remove languages without ISO:
asjp16.dm <- asjp16.dm[ !is.na(lgs.to.keep), !is.na(lgs.to.keep) ]; lgs.to.keep <- lgs.to.keep[ !is.na(lgs.to.keep) ];
# Remove duplicated languages:
asjp16.dm <- asjp16.dm[ !duplicated(lgs.to.keep), !duplicated(lgs.to.keep) ]; lgs.to.keep <- lgs.to.keep[ !duplicated(lgs.to.keep) ];
# Remove languages with the same ISO code:
asjp16.dm <- asjp16.dm[ !duplicated(lgs$ISO[lgs.to.keep]), !duplicated(lgs$ISO[lgs.to.keep]) ]; lgs.to.keep <- lgs.to.keep[ !duplicated(lgs$ISO[lgs.to.keep]) ];
# Set the row and col names to the language ISO:
rownames(asjp16.dm) <- (colnames(asjp16.dm) <- lgs$ISO[lgs.to.keep]);
# Make the matrix a full distances matrix:
asjp16.dm <- asjp16.dm + t(asjp16.dm);
# And normalize it between 0 and 1:
asjp16.dm.min <- min(as.numeric(asjp16.dm)); asjp16.dm.max <- max(as.numeric(asjp16.dm)); asjp16.dm <- (asjp16.dm - asjp16.dm.min)/(asjp16.dm.max - asjp16.dm.min);

if( do.save.as.csv )
{
  # Save this matrix as a CSV file:
  cat( "Writing as CSV (TAB-separated) file)...\n", file="./conversion.log", append=TRUE );
  write.table( asjp16.dm, "./asjp16-dists.csv", quote=FALSE, sep="\t" );
}

if( do.save.as.RData )
{
  # Save as an R object:
  cat( "Writing as R object...\n", file="./conversion.log", append=TRUE );
  save( asjp16.dm, file = "./asjp16-dists.RData",compress="xz" );
}

cat( "Conversion DONE!\n\n", file="./conversion.log", append=TRUE );


# Compare ASJP15 and ASJP16 distances:
load( "./asjp15-dists.RData" ); load( "./asjp16-dists.RData" );
dim(asjp15.dm); dim(asjp16.dm); # 3761 vs 3932 languages
shared.lgs <- sort( intersect( rownames(asjp15.dm), rownames(asjp16.dm) ) ); # 3536 shared ISO codes
asjp15.dm.shared <- asjp15.dm[ shared.lgs, shared.lgs ]; asjp16.dm.shared <- asjp16.dm[ shared.lgs, shared.lgs ]; # the shared part
library(vegan);
mantel( asjp15.dm.shared, asjp16.dm.shared, method="pearson", permutations=100 ); # r=0.98, p=0.0099

