#########################################################################################################################################################
#
# Import, process and export in a standardized Newick format various language classifications (currently WALS, Ethnologue, Glottolog and AUTOTYP)
# and add branch length using a variety or methods (none, constant, proportional, grafen, nnls, ga and nj; see script "TreeHelperFunctions.r" for info)
#
# Copyright (C) 2013-2015  Dan Dediu
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
##########################################################################################################################################################

######################################################################################
#
# Global parametres controlling the various processing, import and export options:
#    tweak as desired (some require lots of computater resources and time)
# More parametres controlling the branch length methods are defined below.
#
######################################################################################
MATCH_CODES       = TRUE; # compute the equivalences between the ISO, WALS, AUTOTYP and GLOTTOLOG codes and save them to file?
PREOPTIMIZE_DISTS = TRUE; # pre-optimize the distance matrices for fast loading when required
COMPUTE_GEO_DISTS = TRUE; # compute the geographic distances between languages?

TRANSFORM_TREES   = TRUE; # transform the trees to the Newick notation (no branch length) and save them?

EXPORT_NEXUS      = TRUE;  # export the trees to a NEXUS file?
EXPORT_NEXUS_TRANSLATE_BLOCK = TRUE; # when exporting NEXUS files, generate a TRANSLATE block (useful when using programs such as BayesTraits that have issues parsing complicated taxa names)?
EXPORT_CSV        = TRUE;  # export the trees to a CSV file?

COMPUTE_BRLEN     = TRUE; # apply the various branch length methods to the Newick no branch length trees?

COMPARE_TREES     = TRUE; # compare (compute the distance between) equivalent trees (within classifications/between branch length methods)?
tree.comparisons  = "../output/tree_comparisons_between_methods.csv"; # file which stores the comparisons

CPU_CORES = 4; # the number of parallel CPU cores to use

quotes="'"; # the quotes to use for the node (language and family) names

#########################
#
#    Tree functions
#
#########################
source( "./FamilyTrees.R" );


######################################################################################
#
# Data needed to map the ISO - WALS - Glottolog codes
#
######################################################################################

if( MATCH_CODES )
{
  # The Ethnologue data contains the ISO and language names:
  ethn.data <- read.table( "../input/ethnologue/LanguageCodes.tab", header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE );
  ethn.data <- ethn.data[ , c("LangID", "Name") ];
  ethn.data$Name <- sapply( as.character(ethn.data$Name), function( s ){ normalize.language.name(s); } );
  
  # The wals data contains the WALS, ISO and Glottolog codes for a subset of languages:
  wals.data <- read.table( "../input/wals/language.csv", header=TRUE, sep=",", quote="\"", stringsAsFactors=FALSE );
  wals.data <- wals.data[ , c("wals_code", "iso_code", "glottocode", "Name", "latitude", "longitude") ];
  wals.data$Name <- sapply( as.character(wals.data$Name), function( s ){ normalize.language.name(s); } );
  
  # The Glottolog data contains the ISO and Glottolog codes for a much larger set of languages:
  library(rjson);
  glott.json <- fromJSON( file="../input/glottolog/glottocodes2iso.json" );
  glott.json <- glott.json$resources; # select the actual resources and not the properties
  glott.data <- lapply( glott.json, function(x){ 
                                        # Look for an ISO code:
                                        isos <- NULL;
                                        if( length(x$identifiers) > 0 )
                                        {
                                          for( i in 1:length(x$identifiers) )
                                          {
                                            if( x$identifiers[[i]]$type == "iso639-3" )
                                            {
                                              isos <- c( isos, x$identifiers[[i]]$identifier );
                                            }
                                          }
                                        }
                                        if( length(isos) > 1 )
                                        {
                                          cat( "Multiple ISOs found for glottocode ", x$id, ": packing them separated by '-'!\n" );
                                          isos <- paste( isos, collapse="-", sep="" );
                                        }
                                        return ( data.frame( "Name"=x$name, 
                                                             "Glottocode"=x$id, 
                                                             "ISO"=ifelse( is.null(isos), NA, isos ), 
                                                             "latitude"=ifelse(is.null(x$latitude),NA,x$latitude),
                                                             "longitude"=ifelse(is.null(x$latitude),NA,x$longitude) ) );
                                      } );
  glott.data <- do.call( rbind, glott.data );
  glott.data$Name <- sapply( as.character(glott.data$Name), function( s ){ normalize.language.name(s); } );
  rm(glott.json);
  
  # The autotyp data contains the ISO, Glottolog and autoty (LID) mappings:
  autotyp.data <- read.table( "../input/autotyp/autotyp-trees.csv", header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE );
  autotyp.data$language <- sapply( as.character(autotyp.data$language), function( s ){ normalize.language.name(s); } );
  
  # The mappings:
  code.mappings <- merge( glott.data, wals.data, by.x=c("Glottocode","ISO"), by.y=c("glottocode","iso_code"), all=TRUE, suffixes=c(".glot",".wals") );
  code.mappings <- merge( code.mappings, autotyp.data, by.x=c("Glottocode","ISO"), by.y=c("glottolog_LID.2014","ISO639.3.2013"), all=TRUE, suffixes=c(".glot",".auto") );
  code.mappings <- merge( code.mappings, ethn.data, by.x=c("ISO"), by.y=c("LangID"), all=TRUE, suffixes=c("",".ethn") );
  # Combine and harmonize the geographic corrdinates:
  # Fill in the missing coordinates in Glottolog from WALS:
  s <- (is.na(code.mappings$latitude.glot) | is.na(code.mappings$longitude.glot));
  code.mappings[s,c("latitude.glot","longitude.glot")] <- code.mappings[s,c("latitude.wals","longitude.wals")];
  # If the Glottolog and WALS coordinates are very different, display the differences and keep the WALS:
  s <- (!is.na(code.mappings$latitude.glot) & !is.na(code.mappings$longitude.glot) & !is.na(code.mappings$latitude.wals) & !is.na(code.mappings$longitude.wals) & 
        (abs(code.mappings$latitude.glot - code.mappings$latitude.wals) >= 1.0 | abs(code.mappings$latitude.glot - code.mappings$latitude.wals) >= 1.0));
  if( sum(s,na.rm=TRUE) > 0 )
  {
    cat( "There are disagreements larger than 1 degree between Glottolog and WALS for ", nrow(unique(code.mappings[s,c("ISO","wals_code","Glottocode")])), " languages (WALS codes): ", paste( unique(code.mappings$wals_code[s]), collapse=",", sep="" ), "; keeping the WALS coordinates!\n" );
    code.mappings[s,c("latitude.glot","longitude.glot")] <- code.mappings[s,c("latitude.wals","longitude.wals")];
  }
  # See how many languages with ISO codes but no geographic coordinates there are:
  code.mappings$ISO[ is.na(code.mappings$ISO) ] <- "";
  s <- ((is.na(code.mappings$latitude.glot) | is.na(code.mappings$longitude.glot)) & code.mappings$ISO != "");
  if( sum(s,na.rm=TRUE) )
  {
    missing.geo <- unique(code.mappings$ISO[s]);
    cat( "There are ", length(missing.geo), " languages with ISO codes but no geographic information: I don't know where to get the coordinates from :(\n" );
  }
  # Data cleaning:
  code.mappings <- code.mappings[ , c( "ISO", "wals_code", "LID", "Glottocode", "Name", "Name.wals", "language", "Name.glot", "latitude.glot", "longitude.glot" ) ]; 
  names(code.mappings) <- c( "ISO", "WALS", "AUTOTYP", "Glottolog", "Name.ethn", "Name.wals", "Name.autotyp", "Name.glottlog", "Latitude", "Longitude" ); 
  code.mappings[is.na(code.mappings)] <- "";
  code.mappings$ISO           <- as.character(code.mappings$ISO);
  code.mappings$WALS          <- as.character(code.mappings$WALS);
  code.mappings$AUTOTYP       <- as.character(code.mappings$AUTOTYP);
  code.mappings$Glottolog     <- as.character(code.mappings$Glottolog);
  code.mappings$Name.ethn     <- as.character(code.mappings$Name.ethn);
  code.mappings$Name.wals     <- as.character(code.mappings$Name.wals);
  code.mappings$Name.autotyp  <- as.character(code.mappings$Name.autotyp);
  code.mappings$Name.glottlog <- as.character(code.mappings$Name.glottlog);
  code.mappings$Latitude      <- as.numeric(code.mappings$Latitude);
  code.mappings$Longitude     <- as.numeric(code.mappings$Longitude);
  code.mappings <- code.mappings[ order( code.mappings$ISO, code.mappings$WALS, code.mappings$AUTOTYP, code.mappings$Glottolog ), ];
  
  # Add a Universal Unique Language Identifier using the conventions introduced here:
  code.mappings <- cbind( code.mappings, "UULID"=paste("[i-",code.mappings$ISO,"][w-",code.mappings$WALS,"][a-",code.mappings$AUTOTYP,"][g-",code.mappings$Glottolog,"]",sep="") );
    
  # Write them to file:
  write.table( code.mappings, "../output/code_mappings_iso_wals_autotyp_glottolog.csv", sep="\t", quote=FALSE, row.names=FALSE );
  code.mappings$AUTOTYP <- as.character(code.mappings$AUTOTYP); code.mappings$AUTOTYP[ is.na(code.mappings$AUTOTYP) ] <- ""; # make sure the AUTOTYP code is treated as character with missing data as ""
} else
{
  code.mappings <- read.table( "../output/code_mappings_iso_wals_autotyp_glottolog.csv", sep="\t", quote="", header=TRUE, stringsAsFactors=FALSE );
  code.mappings$AUTOTYP <- as.character(code.mappings$AUTOTYP); code.mappings$AUTOTYP[ is.na(code.mappings$AUTOTYP) ] <- ""; # make sure the AUTOTYP code is treated as character with missing data as ""
}


######################################################################################
#
# Load (and possibly compute) various distance matrices to be used for computing
# the branch lengths.
#
######################################################################################

# Pre-optimize and cache the distance matrices for on-demand loading only when really needed:
if( PREOPTIMIZE_DISTS )
{
  # ASJP16 (the row and column names are the iso codes):
  load( "../input/distances/ASJP/asjp16-dists.RData" ); # asjp16.dm
  d <- asjp16.dm; save(d, file="../output/preotimized-distances/asjp16-dm.RData", compress="xz"); rm(asjp16.dm,d); # cache it pre-optimized for on-demand loading
  
  # WALS (the row and column names are the wals codes); gowler and manhattan are quite similar; using both NAs and imputation of NAs with the mode (per variable):
  load( "../input/distances/WALS/wals-gower-dm.RData" ); wals.gower.dm <- as.matrix(wals.gower.dm); # wals.gower.dm (with missing data)
  d <- wals.gower.dm; save(d, file="../output/preotimized-distances/wals-gower-dm.RData", compress="xz"); rm(wals.gower.dm,d); # cache it pre-optimized for on-demand loading
  
  load( "../input/distances/WALS/wals-gower-mode-dm.RData" ); wals.gower.mode.dm <- as.matrix(wals.gower.mode.dm); # wals.gower.mode.dm (with mode imputation for missing data)
  d <- wals.gower.mode.dm; save(d, file="../output/preotimized-distances/wals-gower-mode-dm.RData", compress="xz"); rm(wals.gower.mode.dm,d); # cache it pre-optimized for on-demand loading
  
  load( "../input/distances/WALS/wals-euclidean-dm.RData" ); wals.euclidean.dm <- as.matrix(wals.euclidean.dm); # wals.euclidean.dm (with missing data)
  d <- wals.euclidean.dm; save(d, file="../output/preotimized-distances/wals-euclidean-dm.RData", compress="xz"); rm(wals.euclidean.dm,d); # cache it pre-optimized for on-demand loading
  
  load( "../input/distances/WALS/wals-euclidean-mode-dm.RData" ); wals.euclidean.mode.dm <- as.matrix(wals.euclidean.mode.dm); # wals.euclidean.mode.dm (with mode imputation for missing data)
  d <- wals.euclidean.mode.dm; save(d, file="../output/preotimized-distances/wals-euclidean-mode-dm.RData", compress="xz"); rm(wals.euclidean.mode.dm,d); # cache it pre-optimized for on-demand loading
  
  # AUTOTYP (the row and column names are the AUTOTYP codes (LID)):
  load("../input/distances/AUTOTYP/autotyp-dist.RData"); # autotyp.dm (with missing data, using ony features with single values per language, courtesy of Balthasar Bickel)
  d <- autotyp.dm; save(d, file="../output/preotimized-distances/autotyp-dm.RData", compress="xz"); rm(autotyp.dm,d); # cache it pre-optimized for on-demand loading
  
  # Geographic distances (the row and column names are the Glottolog codes):
  if( COMPUTE_GEO_DISTS )
  {
    # Compute the great-circle geographic distances between all languages with geographic coordinates
    library(geosphere);
    # The languages with geographic coordinates:
    s <- (!is.na(code.mappings$Latitude) & !is.na(code.mappings$Longitude)); 
    # of which those that don't have a specific type of code:
    s.ISO       <- s & (code.mappings$ISO == "");
    s.WALS      <- s & (code.mappings$WALS == "");
    s.AUTOTYP   <- s & (code.mappings$AUTOTYP == "");
    s.Glottolog <- s & (code.mappings$Glottolog == "");
    # Decide which code to use so that as many as possible of the languages with geographic info are covered by the coding scheme:
    code.to.use <- c("ISO","WALS","AUTOTYP","Glottolog")[ which.min( c( sum(s.ISO,na.rm=TRUE), sum(s.WALS,na.rm=TRUE), sum(s.AUTOTYP,na.rm=TRUE), sum(s.Glottolog,na.rm=TRUE) ) ) ];
    s.use <- s & !(get( paste( "s.", code.to.use, sep="" ) )); # obtain (s & (code != ""))
    tmp <- code.mappings[s.use,c(code.to.use,"Longitude","Latitude")];
    cat( "For the geographic distances (ellipsoid) the row and column names use the ", code.to.use, " codes.\n" );
    # Some codes (with geographic coordinates) are duplicated: keep the frst entries only:
    tmp <- tmp[!duplicated(tmp[,1]),];
    
    ## Plot these languages:
    #library(maps);
    #X11(10,10);
    #map("world");
    #points( tmp$Longitude, tmp$Latitude, cex=0.5, col=grey(0.5) );
    
    # Compute the distances:
    m <- as.matrix( tmp[,c("Longitude","Latitude")] ); rownames(m) <- tmp[,1];
    system.time( geo.dm <- geosphere::distm( m, fun=distMeeus ) ); # much faster than distVincentyEllipsoid but still a very decent estimate of geographic distance
    rownames(geo.dm) <- colnames(geo.dm) <- rownames(m);
    
    rm(m); rm(tmp);
    
    save( geo.dm, file="../input/distances/geo-great-circle-ellipsoid-dists.RData", compress="xz" );
  } else
  {
    load( file="../input/distances/geo-great-circle-ellipsoid-dists.RData" ); # geo.dm 
    # Normalize it between 0 and 1 by dividing it with the maximum distance (which is very close to the maximum possible distance on Earth of about 21,000km):
    geo.dm <- geo.dm / max(as.numeric(geo.dm),na.rm=TRUE);
    d <- geo.dm; save(d, file="../output/preotimized-distances/geo-dm.RData", compress="xz"); rm(geo.dm,d); # cache it pre-optimized for on-demand loading
  }
  
  # The MG2015 [Maurits, L. & Griffiths, T.L. (2015) PNAS 111 (37):13576--13581] distances
  # there's one distance matrix for each of the four classifications (WALS, Ethnologue, Glottolog and AUTOTYP)
  # and each has a rwo & column names the full language names and code -> we need to replace them by the appropriate codes
  .change.row.col.names.to.codes <- function(m, code=c("iso","wals","autotyp","glottolog"))
  {
    if( is.null(m) || !inherits(m,"matrix") || nrow(m) != ncol(m) || any(rownames(m) != colnames(m)) )
    { 
      stop("Problems with the distance matrix!\n"); 
      return (NULL); 
    } else
    {
      # Change the row and column names to the WALS codes:
      new.names <- vapply(rownames(m), function(s)
        { 
          tmp <- extract.name.and.codes(s); 
          if( is.null(tmp[[code]]) || length(tmp[[code]]) != 1 || tmp[[code]]=="" )
          { 
            warning(paste0("Language ",s," must have a single ", code, " code\n")); 
            return (s); 
          } else 
          {
            return(tmp[[code]]); 
          }
        }, character(1));
      return (new.names);
    }
  }
  
  load("../input/distances/MG2015/MG2015-wals-alpha=0.69.RData"); rownames(MG2015.WALS) <- colnames(MG2015.WALS) <-.change.row.col.names.to.codes(MG2015.WALS, "wals");  # MG2015.WALS
  d <- MG2015.WALS; save(d, file="../output/preotimized-distances/mg2015-wals-dm.RData", compress="xz"); rm(MG2015.WALS,d); # cache it pre-optimized for on-demand loading
  
  load("../input/distances/MG2015/MG2015-ethnologue-alpha=0.69.RData"); rownames(MG2015.Ethnologue) <- colnames(MG2015.Ethnologue) <-.change.row.col.names.to.codes(MG2015.Ethnologue, "iso");  # MG2015.Ethnologue
  d <- MG2015.Ethnologue; save(d, file="../output/preotimized-distances/mg2015-ethnologue-dm.RData", compress="xz"); rm(MG2015.Ethnologue,d); # cache it pre-optimized for on-demand loading
  
  load("../input/distances/MG2015/MG2015-glottolog-alpha=0.69.RData"); rownames(MG2015.Glottolog) <- colnames(MG2015.Glottolog) <-.change.row.col.names.to.codes(MG2015.Glottolog, "glottolog");  # MG2015.Glottolog
  d <- MG2015.Glottolog; save(d, file="../output/preotimized-distances/mg2015-glottolog-dm.RData", compress="xz"); rm(MG2015.Glottolog,d); # cache it pre-optimized for on-demand loading
  
  load("../input/distances/MG2015/MG2015-autotyp-alpha=0.69.RData"); rownames(MG2015.AUTOTYP) <- colnames(MG2015.AUTOTYP) <-.change.row.col.names.to.codes(MG2015.AUTOTYP, "autotyp");  # MG2015.AUTOTYP
  d <- MG2015.AUTOTYP; save(d, file="../output/preotimized-distances/mg2015-autotyp-dm.RData", compress="xz"); rm(MG2015.AUTOTYP,d); # cache it pre-optimized for on-demand loading
  
  gc(); # call the garbage collector to make sure the space is freeded after this memory-hungry step
}


######################################################################################
#
# Global parametres controlling the branch length methods to apply
#
######################################################################################

# the branch length methods to be used; can be "all", "none", or any subset of {"constant", "proportional", "grafen", "nnls", "ga", "nj"}
CLASSIFICATIONS = c("wals", "ethnologue", "glottolog", "autotyp");
METHODS   = c("constant", "proportional", "grafen", "nnls", "nj", "ga"); 
CONSTANT  = 1.0; # the positive constant required by some methods
DISTS.CODES = read.table(text="Distance              ShortName  File                                                          Code      \n
                               asjp16                asjp       ../output/preotimized-distances/asjp16-dm.RData               iso       \n
                               wals(gower)           w:g        ../output/preotimized-distances/wals-gower-dm.RData           wals      \n
                               wals(gower,mode)      w:gm       ../output/preotimized-distances/wals-gower-mode-dm.RData      wals      \n
                               wals(euclidean)       w:e        ../output/preotimized-distances/wals-euclidean-dm.RData       wals      \n
                               wals(euclidean,mode)  w:em       ../output/preotimized-distances/wals-euclidean-mode-dm.RData  wals      \n
                               autotyp               atyp       ../output/preotimized-distances/autotyp-dm.RData              autotyp   \n
                               mg2015(wals)          mg:w       ../output/preotimized-distances/mg2015-wals-dm.RData          wals      \n
                               mg2015(ethnologue)    mg:e       ../output/preotimized-distances/mg2015-ethnologue-dm.RData    iso       \n
                               mg2015(glottolog)     mg:g       ../output/preotimized-distances/mg2015-glottolog-dm.RData     glottolog \n
                               mg2015(autotyp)       mg:a       ../output/preotimized-distances/mg2015-autotyp-dm.RData       autotyp   \n
                               geo                   geo        ../output/preotimized-distances/geo-dm.RData                  glottolog \n",
                         sep="", quote="", header=TRUE, stringsAsFactors=FALSE); # the distances required by some methods
              

######################
#
# WALS trees
#
######################

if( TRANSFORM_TREES )
{
  # WALS languages info:
  wals.data <- read.table( "../input/wals/language.csv", header=TRUE, sep=",", quote="\"", stringsAsFactors=FALSE );
  wals.data <- wals.data[, c("wals_code", "iso_code", "glottocode", "Name", "latitude", "longitude", "genus", "family" ) ];
  
  # Replace troublesome characters in language and family names:
  wals.data$Name <- sapply( as.character(wals.data$Name), function( s ){ normalize.language.name(s); } );
  wals.data$family <- sapply( as.character(wals.data$family), function( s ){ normalize.language.name(s); } );
  wals.data$genus <- sapply( as.character(wals.data$genus), function( s ){ normalize.language.name(s); } );
  
  # Create the WALS trees:
  roots <- NULL;
  for( i in 1:nrow(wals.data) )
  {
    # Process this entry:
    if( !is.na(wals.data$family[i]) && !(wals.data$family[i] %in% c("other","Unclassified","Creole","Mixed language")) )
    {
      autotypcode <- unique( as.character( code.mappings$AUTOTYP[ code.mappings$WALS == as.character(wals.data$wals_code[i]) ] ) );
      if( length(autotypcode) > 1 )
      {
        # Get rid of empty codes:
        autotypcode <- autotypcode[ autotypcode != "" ];
        if( length(autotypcode) == 0 )
        {
          autotypcode = ""; # all were ""
        } else
        {
          autotypcode <- paste( autotypcode, collapse="-", sep="" ); # in principle it can happen to have more than one AUTOTYP code
        }
      }
      roots <- add.language.to.languageclassification( roots, 
                                                       path=c(create.name.and.codes(wals.data$family[i], iso="", wals="", autotyp="", glottocode="", quotes=quotes),  
                                                              create.name.and.codes(wals.data$genus[i],  iso="", wals="", autotyp="", glottocode="", quotes=quotes),
                                                              create.name.and.codes(wals.data$Name[i], 
                                                                                    iso=       as.character(wals.data$iso_code[i]), 
                                                                                    wals=      as.character(wals.data$wals_code[i]), 
                                                                                    autotyp=   autotypcode, 
                                                                                    glottocode=as.character(wals.data$glottocode[i]), 
                                                                                    quotes) ),
                                                       brlens=NULL );
    }
  }
  attr(roots, "classif.name") <- "wals";
  # Fix the family and language names (get rid of the codes):
  for( i in 1:length(roots) )
  {
    roots[[i]] <- .fix.names(roots[[i]], quotes=quotes);
    attr(roots[[i]], "tree.name") <- str_trim(strsplit(strsplit(get.name(roots[[i]]),quotes,fixed=TRUE)[[1]][2],"[",fixed=TRUE)[[1]][1]);
  }
  # Save them to file:
  export.languageclassification( roots, dir.name="../output", classification="wals", export.nexus=EXPORT_NEXUS, nexus.translate.block=EXPORT_NEXUS_TRANSLATE_BLOCK, export.csv=EXPORT_CSV );
  
  rm(wals.data);
}


######################
#
# Ethnologue trees
#
######################

if( TRANSFORM_TREES )
{
  # The Enthologue trees: we need to download and parse them from the Ethnologue website:
  fams.file <- paste("../input/ethnologue/fams.html",sep="");
  fams.URL <- "http://www.ethnologue.com/browse/families";
  if( !file.exists(fams.file) || file.info(fams.file)$size < 1 )
  {
    options(timeout=180); # Ethnologue is a slow site for some reason...
    options(HTTPUserAgent="Mozilla/5.0 (X11; Linux x86_64; rv:17.0) Gecko/20130917 Firefox/17.0"); # emulate Firefox
    if( download.file( fams.URL, fams.file, quiet=T ) != 0 )
    {
      stop( "Error retrieving the Ethnologue main families page !\n" );
    }
  }
  # Load the file:
  HTML.content <- readLines( fams.file, warn=FALSE ); # read the file:
  # Find the links to the families themselves:
  fams.links <- HTML.content[ grep( "/subgroups/", HTML.content, fixed=TRUE ) ];
  # Extract the links:
  fams.links <- sapply( fams.links, function(s){ 
                                                 tmp <- strsplit( s, "<a href=\"/subgroups/", fixed=TRUE )[[1]];
                                                 if( length(tmp) < 2 )
                                                 {
                                                   cat( "Error retrieving the URL for ", s, "\n" );
                                                   return (NULL);
                                                 }
                                                 tmp <- strsplit( tmp[2], "\">", fixed=TRUE )[[1]];
                                                 if( length(tmp) < 1 )
                                                 {
                                                   cat( "Error retrieving the URL for ", s, "\n" );
                                                   return (NULL);
                                                 }
                                                 return (tmp[1]);
                                               } );
  
  # Load and process one language family
  fetch.ethnologue.family <- function( family )
  {
    if( is.na( family ) )
    {
      return( NULL );
    }
    else
    {
      # The base URL for families:
      fam.base.URL <- "http://www.ethnologue.com/subgroups/";
      fam.URL <- paste( fam.base.URL, family, sep="" );
      fam.file <- paste("../input/ethnologue/fam-",family,".html",sep="");
      if( !file.exists(fam.file) || file.info(fam.file)$size < 1 )
      {
        options(timeout=180); # Ethnologue is a slow site for some reason...
        options(HTTPUserAgent="Mozilla/5.0 (X11; Linux x86_64; rv:17.0) Gecko/20130917 Firefox/17.0"); # emulate Firefox
        if( download.file( fam.URL, fam.file, quiet=T ) != 0 )
        {
          cat( "Error retrieving family ", family, " !\n" );
          return( NULL );
        }
      }
      
      # Use the downloaded file:
      HTML.content <- readLines( fam.file, warn=F ); # read the file:
      
      # Get the family official name:
      tag.string <- "<h1 class=\"title\" id=\"page-title\">";
      matches <- grep( tag.string, HTML.content, fixed=TRUE ); 
      if( length(matches) != 1 )
      {
        cat( "Error locating family name for ", family, "\n" )
        return ( NULL );
      } else
      {
        # Get the family name:
        family.name <- strsplit( strsplit( HTML.content[ matches ], tag.string, fixed=TRUE )[[1]][2], "</h1>", fixed=TRUE )[[1]][1];
      }
      
      # The tree in preorder with special symbols ( and ) representing an increse/decrease in level
      # the names are "node" for an internal node, "leaf" for a leaf, and "special" for ( and ):
      tree <- c("node"=create.name.and.codes( family.name, iso="", wals="", autotyp="", glottocode="", quotes=quotes));
      
      # Preprocess the HTML text to make sure each of the important tags appear only once per line (i.e., insert newlines after each </a> and </ul>):
      HTML.content.new <- gsub("</a>", "</a>\n", HTML.content[(matches+1):length(HTML.content)], fixed=TRUE);
      HTML.content.new <- gsub("</ul>", "</ul>\n", HTML.content.new, fixed=TRUE);
      HTML.content.new <- unlist(strsplit(HTML.content.new, "\n", fixed=TRUE));
      
      # Build the tree by adding the subfamilies and languages as they appear in the file:
      for( i in 1:length(HTML.content.new) )
      {
        if( length( grep( "<a href=\"/subgroups/", HTML.content.new[i], fixed=TRUE ) ) == 1 )
        {
          # This seems to be a subgroup: extract its name:
          subgroup <- strsplit( strsplit( strsplit( HTML.content.new[i], "<a href=\"/subgroups/", fixed=TRUE )[[1]][2], "\">", fixed=TRUE )[[1]][2], "</a>", fixed=TRUE )[[1]][1];
          subgroup <- str_trim( strsplit( subgroup, "\\([[:digit:]]+\\)" )[[1]][1] );
          if( is.null(subgroup) || subgroup == "" )
          {
            cat( "Error extracting subgroup from ", HTML.content.new[i], "\n" );
            return ( NULL );
          }
          # Add this to the tree:
          tree <- c(tree, "node"=create.name.and.codes(subgroup, iso="", wals="", autotyp="", glottocode="", quotes=quotes)); 
        } else if( length( grep( "<span class=\"field-content\">", HTML.content.new[i], fixed=TRUE ) ) == 1 && length( grep( "<a href=\"/language/", HTML.content.new[i], fixed=TRUE ) ) == 1 )
        {
          # Seems to be language: extract it and its iso:
          lg.name <- strsplit( HTML.content.new[i], "<span class=\"field-content\">", fixed=TRUE )[[1]][2];
          lg.name <- str_trim( strsplit( lg.name, "<a href=\"/language/", fixed=TRUE )[[1]][1] );
          lg.iso <- strsplit( HTML.content.new[i], "<a href=\"/language/", fixed=TRUE )[[1]][2];
          lg.iso <- strsplit( lg.iso, "\">[", fixed=TRUE )[[1]][2];
          lg.iso <- strsplit( lg.iso, "]</a>", fixed=TRUE )[[1]][1];
          if( is.null(lg.name) || lg.name == "" ||  is.null(lg.iso) || lg.iso == "" )
          {
            cat( "Error extracting language from ", HTML.content.new[i], "\n" );
            return ( NULL );
          }
          # Get the other codes as well:
          wals <- autotyp <- glottocode <- "";
          s <- code.mappings[ code.mappings$ISO == lg.iso, ];
          if( !is.null(s) && nrow(s) > 0 )
          {
            wals       <- unique(s$WALS);
            autotyp    <- unique(s$AUTOTYP);
            glottocode <- unique(s$Glottolog);
          }
          # Add this language to the tree:
          tree <- c(tree, "leaf"=create.name.and.codes(lg.name, iso=lg.iso, wals=wals, autotyp=autotyp, glottocode=glottocode, quotes=quotes));
        } else if( length( grep( "<ul>", HTML.content.new[i], fixed=TRUE ) ) == 1 )
        {
          # Subgroup starts:
          tree <- c(tree, "special"="(");
        } else if( length( grep( "</ul>", HTML.content.new[i], fixed=TRUE ) ) == 1 )
        {
          # Subgroup ends:
          tree <- c(tree, "special"=")");
        } else if( length( grep( "</section>", HTML.content.new[i], fixed=TRUE ) ) == 1 )
        {
          # The end of the tree section
          break;
        }
      }
      
      # Build the tree reusing as much code as possible: simply add the full path for each language:
      root <- NULL;
      cur.path <- tree[1];
      for( i in 2:length(tree) )
      {
        #cat(paste0("Tree[",i,"]=",tree[i],": "));
        if( names(tree)[i] == "special" && tree[i] == "(" )
        {
          # Simply ignore, we're going up one level
        } else if( names(tree)[i] == "special" && tree[i] == ")" )
        {
          if( i < length(tree) && names(tree)[i+1] == "special" && tree[i+1] == "(" ) next; # artifact of the way languages and sub-families are mixed in the HTML file
          # Down one level, discard the last node in the path:
          cur.path <- cur.path[ -length(cur.path) ];
          #cat("cur.path=",paste0(cur.path,collapse=", "));
        } else if( names(tree)[i] == "leaf" )
        {
          # Add the full path to the tree:
          root <- add.tree.path(root, c(cur.path, tree[i]));
          #cat(paste0("add.tree.path(",paste0(c(cur.path, tree[i]),collapse=", "),")"));
        } else if( names(tree)[i] == "node" )
        {
          # Add the node to the path:
          cur.path <- c(cur.path, tree[i]);
          #cat("cur.path=",paste0(cur.path,collapse=", "));
        }
        #cat("\n");
      }
      attr(root, "tree.name") <- family.name;
      
      return( root );
    }
  }
  fetch.ethnologue.family <- cmpfun(fetch.ethnologue.family);
  
  # Fetch the families:
  roots <- list();
  for( i in 1:length(fams.links) )
  {
    cat( "Processing Ethnologue family ", fams.links[i], "...\n" );
    root <- fetch.ethnologue.family( fams.links[i] );
    if( is.null(root) )
    {
      cat( "  >>>  Error fetching Ethnologue family ", fams.links[i], "...\n" );
      next;
    }
    roots[[ length(roots) + 1 ]] <- .fix.names( root, quotes=quotes );
  }
  # Save them to file:
  attr(roots, "classif.name") <- "Ethnologue";
  class(roots) <- append("languageclassification", class(roots));
  export.languageclassification( roots, dir.name="../output", classification="ethnologue", export.nexus=EXPORT_NEXUS, nexus.translate.block=EXPORT_NEXUS_TRANSLATE_BLOCK, export.csv=EXPORT_CSV );
  
  rm(HTML.content); rm(fams.links);
}


######################
#
# Glottolog trees
#
######################

if( TRANSFORM_TREES )
{
  # The Glottolog trees: parse them and export them in the Family\tTree format, also updating the convention on the language codes:
  glott.fams <- readLines( "../input/glottolog/tree-glottolog-newick.txt" );
  mapping <- code.mappings[,c("ISO","WALS","AUTOTYP","Glottolog")]
  roots <- NULL;
  for( i in 1:length(glott.fams) )
  {
    # Pre-process the tree by replacing any "\\'" within langauge names by "`", and taking care of names with ":" by replacing it with "-":
    gfam <- gsub( "\\'", "`", str_trim(glott.fams[i]), fixed=TRUE );
    tmp <- strsplit( gfam, "'", fixed=TRUE )[[1]];
    k <- grep( "[^:]+:", tmp );
    if( length(k) > 0 )
    {
      # ":" in names!
      for( j in k ) tmp[j] <- gsub( ":", "-", tmp[j], fixed=TRUE );
      gfam <- paste( tmp, collapse="'", sep="" );
    }
    
    # Make sure the Newick is enclosed between '(' and ');':
    if( substr(gfam,1,2) != "(" || substr(gfam,nchar(gfam)-1,nchar(gfam)) != ");" )
    {
      if( substr(gfam,nchar(gfam),nchar(gfam)) == ";" ) gfam <- paste0("(",substr(gfam,1,nchar(gfam)-1),");") else gfam <- paste0("(",gfam,");");
    }
    
    # Parse the tree:
    root <- familytree(gfam);
    # Set the family name (the root's only child):
    attr(root, "tree.name") <- str_trim(strsplit(strsplit(.get.node.name(root,root$edge[root$edge[,1] == find.root(root), 2][1]),"[",fixed=TRUE)[[1]][1],quotes,fixed=TRUE)[[1]][2]);
    # Convert the Glottolog convention to our conventions in what concerns the codes:
    root <- .convert.glottolog.convention(root, mapping);
    # Erase the branch length info:
    root <- compute.brlen(root, "none")$tree;
    
    # Fix the names and add it to the collection:
    roots[[ length(roots) + 1 ]] <- .fix.names(root, quotes);
  }
  # Save them to file:
  attr(roots, "classif.name") <- "Glottolog";
  class(roots) <- append("languageclassification", class(roots));
  export.languageclassification( roots, dir.name="../output", classification="glottolog", export.nexus=EXPORT_NEXUS, nexus.translate.block=EXPORT_NEXUS_TRANSLATE_BLOCK, export.csv=EXPORT_CSV );
  
  rm(glott.fams); rm(mapping);
}


######################
#
# AUTOTYP trees
#
######################

if( TRANSFORM_TREES )
{
  # AUTOTYP data:
  autotyp.data <- read.table( "../input/autotyp/autotyp-trees.csv", header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE );
  
  # Replace troublesome characters in language, and branch names:
  autotyp.data$language <- sapply( as.character(autotyp.data$language), function( s ){ if( !is.na(s) ) { normalize.language.name(s); } else { NA; } } );
  autotyp.data$lsbranch <- sapply( as.character(autotyp.data$lsbranch), function( s ){ if( !is.na(s) ) { normalize.language.name(s); } else { NA; } } );
  autotyp.data$ssbranch <- sapply( as.character(autotyp.data$ssbranch), function( s ){ if( !is.na(s) ) { normalize.language.name(s); } else { NA; } } );
  autotyp.data$sbranch  <- sapply( as.character(autotyp.data$sbranch),  function( s ){ if( !is.na(s) ) { normalize.language.name(s); } else { NA; } } );
  autotyp.data$mbranch  <- sapply( as.character(autotyp.data$mbranch),  function( s ){ if( !is.na(s) ) { normalize.language.name(s); } else { NA; } } );
  autotyp.data$stock    <- sapply( as.character(autotyp.data$stock),    function( s ){ if( !is.na(s) ) { normalize.language.name(s); } else { NA; } } );
  
  # Replace the NAs in Glottolog and ISO codes by "":
  autotyp.data$glottolog_LID.2014[ is.na(autotyp.data$glottolog_LID.2014) ] <- "";
  autotyp.data$ISO639.3.2013[ is.na(autotyp.data$ISO639.3.2013) ] <- "";
  
  # Create the AUTOTYP trees:
  roots <- NULL;
  for( i in 1:nrow(autotyp.data) )
  {
    # Process this entry:
    if( !is.na(autotyp.data$stock[i]) )
    {
      # Build the path; the stock is always given:
      path <- create.name.and.codes(autotyp.data$stock[i], iso="", wals="", autotyp="", glottocode="", quotes=quotes); 
      # For the mbranch, sbranch, ssbranch & lsbranch things are a bit more complex in the sense that sometimes an intermediate level is skipped (NA) but deeper levels are given: add them with later tree-wide disambiguation:
      defined.levels <- !is.na(autotyp.data[ i, c( "mbranch", "sbranch", "ssbranch", "lsbranch" ) ]);
      if( any(defined.levels) )
      {
        # There is such a level defined
        deepest.defined.level <- max(which( !is.na(autotyp.data[ i, c( "mbranch", "sbranch", "ssbranch", "lsbranch" ) ]) ));
        mbranch.column <- which(names(autotyp.data)=="mbranch")-1;
        for( j in 1:deepest.defined.level )
        {
          if( is.na(autotyp.data[ i, j+mbranch.column ]) )
          {
            path <- c(path, create.name.and.codes(paste0("Group",j), iso="", wals="", autotyp="", glottocode="", quotes=quotes));
          } else
          {
            path <- c(path, create.name.and.codes(autotyp.data[ i, j+mbranch.column ], iso="", wals="", autotyp="", glottocode="", quotes=quotes));
          }
        }
      }
      # Add the language, also finding the other matching code(s):
      s <- which(code.mappings$AUTOTYP == as.character(autotyp.data$LID[i]));
      path <- c(path, create.name.and.codes(autotyp.data$language[i], 
                                            iso=as.character(autotyp.data$ISO639.3.2013[i]), 
                                            wals=ifelse( length(s) == 0, "", paste(code.mappings$WALS[s],collapse="-",sep="") ), 
                                            autotyp=as.character(autotyp.data$LID[i]),
                                            glottocode=as.character(autotyp.data$glottolog_LID.2014[i])));
      
      roots <- add.language.to.languageclassification( roots, 
                                                       path=path,
                                                       brlens=NULL );
    }
  }
  attr(roots, "classif.name") <- "autotyp";
  # Fix the family and language names (get rid of the codes):
  for( i in 1:length(roots) )
  {
    roots[[i]] <- .fix.names(roots[[i]], quotes=quotes);
    attr(roots[[i]], "tree.name") <- str_trim(strsplit(strsplit(get.name(roots[[i]]),quotes,fixed=TRUE)[[1]][2],"[",fixed=TRUE)[[1]][1]);
  }
  # Save them to file:
  export.languageclassification( roots, dir.name="../output", classification="autotyp", export.nexus=EXPORT_NEXUS, nexus.translate.block=EXPORT_NEXUS_TRANSLATE_BLOCK, export.csv=EXPORT_CSV );
  
  rm(autotyp.data);
}


#################################
#
# Compute the branch lengths
#
#################################

if( COMPUTE_BRLEN )
{
  # Apply various methods for computing branch length to a classification:
  compute.brlength.trees <- function( dir.name,                   # the directory where the original trees are and where the resulting trees will be saved
                                      classification,             # the classification: there must be a trees file in ./output/classification/classification-newick.csv
                                      export.nexus=TRUE,          # export the resulting trees to a NEXUS file?
                                      nexus.translate.block=TRUE, # when exporting to NEXUS, also generate a translate block?
                                      export.csv=TRUE,            # export the resulting trees to a CSV file?
                                      methods="all",              # the branch length methods to be used; can be "all", "none", or any subset of {"constant", "proportional", "grafen", "nnls", "ga", "nj"}
                                      constant=1.0,               # the positive constant required by the methods "constant" and "proportional"
                                      dists.codes=NULL,           # the data.frame of distance matrices to be used by "nnls", "ga" and "nj", and the codes used as the distance matrices' row and column names
                                      replace.NA.brlen.with=NA,   # some methods produce NA branch length (if the actual brlen cannot be estimated from the data): replace it by another numeric value?
                                      restore.collapsed.singles=TRUE, # some standard methods need to collapse singles, but we can restore them
                                      parallel.mc.cores=4,
                                      quotes="'"
  )
  {
    # Helper function: given a distances matrix with row and column names as language codes and which code to use, generate the corresponding submatrix with row and column full language names:
    dists.submatrix <- function( roots, dist.mat, code=c("iso","wals","autotyp","glottolog") )
    {
      # Sanity checks:
      if( is.null(roots) || length(roots) == 0 || is.null(dist.mat) || nrow(dist.mat) == 0 || ncol(dist.mat) == 0 )
      {
        return (NULL);
      }
      
      # Adapt the distance methods for these trees (keep only the languages in the intersection and generate the full names):
      lgs <- data.frame( "full"=unlist( lapply( roots, function(root){ extract.languages( root ); } ) ), "name"=NA, "iso"=NA, "wals"=NA, "autotyp"=NA, "glottolog"=NA, stringsAsFactors=FALSE );
      for( i in 1:nrow(lgs) )
      {
        tmp <- extract.name.and.codes( lgs$full[i], quotes );
        lgs$name[i]      <- normalize.language.name(tmp$name);
        lgs$iso[i]       <- paste0(tmp$iso,collapse="-");
        lgs$wals[i]      <- paste0(tmp$wals,collapse="-");
        lgs$autotyp[i]   <- paste0(tmp$autotyp,collapse="-");
        lgs$glottolog[i] <- paste0(tmp$glottolog,collapse="-");
      }
      
      # Expand the required code in the long format, keeping only the actually defined ones:
      codes <- lapply( 1:nrow(lgs), function(i){ tmp <- ifelse( lgs[i,code] == "", "", strsplit(lgs[i,code],"-",fixed=TRUE)[[1]] ); data.frame( "name"=lgs$full[i], "code"=tmp ); } );
      codes <- do.call( rbind, codes ); 
      codes <- codes[ codes$code != "", ];
      
      # Overlap santiy check:
      if( length( intersect( codes$code, rownames(dist.mat) ) ) == 0 )
      {
        return (NULL);
      }
      
      # The shared codes:
      shared.codes <- codes[ codes$code %in% rownames(dist.mat), ];
      if( is.null(shared.codes) || length(shared.codes) == 0 )
      {
        return (NULL);
      }
      # Select the submatrix corresponding to these codes
      dm <- dist.mat[ as.character(shared.codes$code), as.character(shared.codes$code), drop=FALSE ];
      if( is.null(dm) || !is.matrix(dm) || nrow(dm) == 0 || ncol(dm) == 0 )
      {
        return (NULL);
      }
      shared.codes <- shared.codes[ shared.codes$code %in% rownames(dm), , drop=FALSE ];
      if( is.null(shared.codes) || nrow(shared.codes) == 0 )
      {
        return (NULL);
      }
      full.language.names <- as.character(shared.codes$name);
      
      # For duplicated languages, keep the code that minimized the matrix-wide amount of missing data:
      duplgs <- duplicated( shared.codes$name ); 
      if( sum(duplgs,na.rm=TRUE) > 0 )
      {
        # The duplicated language names:
        duplgnames <- unique(shared.codes$name[ duplgs ]);
        # The tried combinations, each composed of an iso code per duplicated language and the resulting number of missing data in the distances matrix:
        lgcombs <- list();
        for( i in 1:length(duplgnames) )
        {
          lgcombs[[i]] <- shared.codes$code[ shared.codes$name==duplgnames[i] ];
        }
        lgcombs <- expand.grid( lgcombs );
        names(lgcombs) <- sapply( 1:length(duplgnames), function(i){ duplgnames[i] } );
        lgcombs <- cbind( lgcombs, "NoNAs"=NA );
        # The non-duplicated language codes:
        nondupcodes <- as.character(shared.codes$code[ !(shared.codes$name %in% duplgnames) ]);
        # Compute the number of NAs for each combination:
        for( i in 1:nrow(lgcombs) )
        {
          # The distances matrix:
          combcodes <- c( nondupcodes, as.character(unlist(lgcombs[i,1:ncol(lgcombs)-1])) );
          d <- dm[ combcodes, combcodes ];
          lgcombs$NoNAs[i] <- sum( is.na(d) );
        }
        if( sum( !is.na(lgcombs$NoNAs) ) > 0 )
        {
          # Pick one that minimizes the number of NAs:
          combcodes <- c( nondupcodes, as.character(unlist(lgcombs[which.min(lgcombs$NoNAs),1:ncol(lgcombs)-1])) );
          dm <- dm[ combcodes, combcodes ];
          # And the corresponding language names:
          full.language.names <- c( as.character(shared.codes$name[ !(shared.codes$name %in% duplgnames) ]), as.character(duplgnames) );
        }
      } 
      # Make the row and column names the full language names:
      rownames(dm) <- colnames(dm) <- full.language.names;
      
      # Return this submatrix:
      return (dm);
    }
    
    # Sanity checks:
    if( constant <= 0 )
    {
      cat( "The constant must be positive\n" );
      return (FALSE);
    }
    
    # Expand the methods to be used:
    if( length(methods) == 1 && methods == "none" )
    {
      # Nothing to do:
      return (TRUE);
    } else if( length(methods) == 1 && methods == "all" )
    {
      # All methods
      methods = c("constant", "proportional", "grafen", "nj", "nnls", "ga");
    } else if( !all( methods %in% c("constant", "proportional", "grafen", "nj", "nnls", "ga"), na.rm=TRUE ) )
    {
      cat( "Unknown methods\n" );
      return (FALSE);
    }   
    
    # Read the roots from file:
    cat( paste( "Reading the trees for classification '", classification, "'...\n", sep="") );
    roots <- languageclassification( classif.name=classification, csv.file=paste(dir.name,"/",classification,"/",classification,"-newick.csv",sep=""), csv.name.column="Family", csv.tree.column="Tree" );
    if( is.null(roots) || length(roots) == 0 || !inherits(roots,"languageclassification") )
    {
      return (FALSE);
    }
    
    # Compute the branch lengths and export them as separate files, one per method (by distance or constant):
    for( method in methods )
    {
      # Decide which constant and/or distances to use depending on the method (and store that in the d.c.k data.frame):
      if( method == "grafen" )
      {
        d.c.k <- data.frame("Distance"=NA, "ShortName"=NA, "File"=NA, "Code"=NA, "k"=NA);
      } else if( method %in% c("constant", "proportional") )
      {
        d.c.k <- data.frame("Distance"=NA, "ShortName"=NA, "File"=NA, "Code"=NA, "k"=constant);
      } else if( method %in% c("nnls", "ga", "nj") )
      {
        d.c.k <- cbind(dists.codes, "k"=NA);
      } else
      {
        # Unknown method; how did I get here?
        next;
      }
      
      # Apply the method:
      for( i in 1:nrow(d.c.k) )
      {
        # Print progress message:
        cat( paste( "\t\tProcesssing the '", classification, "' trees with method '", as.character(method), "'", 
                    ifelse( is.na(d.c.k$k[i]), "", paste( " with constant=", sprintf("%.02f",d.c.k$k[i]), sep="" ) ), 
                    ifelse( is.na(d.c.k$Distance[i]), "", paste( " with distance matrix '", d.c.k$Distance[i], "'", sep="" ) ),
                    "...\n", sep="" ) );
        
        # Load and adapt the distance matrix:
        if( !is.na(d.c.k$Distance[i]) )
        {
          load(d.c.k$File[i]); if( is.null(d) || class(d) != "matrix" || nrow(d) != ncol(d) || nrow(d) < 1 ){ cat("Distance ", d.c.k$Distance[i], " is malformed\n"); return (FALSE); }
          d <- dists.submatrix(roots, d,  d.c.k$Code[i]); 
          if( is.null(d) || class(d) != "matrix" || nrow(d) != ncol(d) || nrow(d) < 1 ){ cat("Distance ", d.c.k$Distance[i], " is empty for classification ", classification, "\n"); return (FALSE); }
        } else
        {
          d <- NULL;
        }
        
        # Apply the method:
        export.languageclassification( roots,
                                       dir.name=dir.name,
                                       classification=classification, 
                                       export.nexus=export.nexus,
                                       nexus.translate.block=nexus.translate.block,
                                       export.csv=export.csv,
                                       method=as.character(method), 
                                       constant=d.c.k$k[i], 
                                       distmatrix=d,
                                       replace.NA.brlen.with=replace.NA.brlen.with, 
                                       restore.collapsed.singles=restore.collapsed.singles,
                                       filename.postfix=paste( "-", method, 
                                                               ifelse( is.na(d.c.k$k[i]), "", paste( "=", sprintf("%.02f",d.c.k$k[i]), sep="") ), 
                                                               ifelse( is.na(d.c.k$Distance[i]), "", paste( "+", d.c.k$Distance[i], sep="") ), 
                                                               sep="" ), 
                                       quotes=quotes, 
                                       parallel.mc.cores=parallel.mc.cores );
        
        # Free up the space:
        rm(d);
      }
    }
    
    return (TRUE);
  }
  compute.brlength.trees <- cmpfun(compute.brlength.trees);

  # Compute the branch lengths for each classification:
  for( classification in CLASSIFICATIONS )
  {
    # Apply only the relevant mg2015 distances:
    dists.codes <- DISTS.CODES[ -grep("mg2015", DISTS.CODES$Distance, fixed=TRUE), ]; dists.codes <- rbind(dists.codes, DISTS.CODES[ DISTS.CODES$Distance == paste0("mg2015(",classification,")"), ]);
    # Compute the branch lengths:
    compute.brlength.trees( dir.name="../output", classification=classification, 
                            export.nexus=EXPORT_NEXUS, nexus.translate.block=EXPORT_NEXUS_TRANSLATE_BLOCK, export.csv=EXPORT_CSV,
                            methods=METHODS, constant=CONSTANT, dists.codes=dists.codes, replace.NA.brlen.with=NA, restore.collapsed.singles=TRUE,
                            parallel.mc.cores=CPU_CORES, quotes=quotes );
  }
}


##################################
#
# Compare trees across methods
#
##################################

if( COMPARE_TREES )
{
  # For a given classification compare all corresponding trees for various methods
  compare.brlength.trees <- function( dir.name,          # the directory where the original trees are and where the resulting trees will be saved
                                      classification,    # the classification: there must be a trees file in ./output/classification/classification-newick.csv
                                      methods="all",     # the branch length methods to be used; can be "all", "none", or any subset of {"constant", "proportional", "grafen", "nnls", "ga", "nj"}
                                      constant=1.0,      # the positive constant required by the methods "constant" and "proportional"
                                      dists.codes=NULL,  # the data.frame of distance matrices to be used by "nnls", "ga" and "nj", and the codes used as the distance matrices' row and column names
                                      parallel.mc.cores=4,
                                      quotes="'"
  )
  {
    # Sanity checks:
    if( constant <= 0 )
    {
      cat( "The constant must be positive\n" );
      return (FALSE);
    }
    
    # Expand the methods to be used:
    if( length(methods) == 1 && methods == "none" )
    {
      # Nothing to do:
      return (TRUE);
    } else if( length(methods) == 1 && methods == "all" )
    {
      # All methods
      methods = c("constant", "proportional", "grafen", "nnls", "ga", "nj");
    } else if( !all( methods %in% c("constant", "proportional", "grafen", "nnls", "ga", "nj"), na.rm=TRUE ) )
    {
      cat( "Unknown methods\n" );
      return (FALSE);
    }  
    
    # The return value is a data.frame which for each family tree and two methods gives the distance(s):
    all.methods <- lapply( methods, function(method)
    {
      # Decide which constant and/or distances to use depending on the method (and store that in the d.c.k data.frame):
      if( method == "grafen" )
      {
        d.c.k <- data.frame("Distance"=NA, "ShortName"=NA, "File"=NA, "Code"=NA, "k"=NA);
      } else if( method %in% c("constant", "proportional") )
      {
        d.c.k <- data.frame("Distance"=NA, "ShortName"=NA, "File"=NA, "Code"=NA, "k"=constant);
      } else if( method %in% c("nnls", "ga", "nj") )
      {
        d.c.k <- cbind(dists.codes, "k"=NA);
      } else
      {
        # Unknown method; how did I get here?
        return (NULL);
      }
      
      ret.val <- lapply( 1:nrow(d.c.k), function(i){
        return(data.frame("Classification"=classification,
                          "Method"=as.character(method),
                          "Constant"=d.c.k$k[i],
                          "Distance"=d.c.k$Distance[i],
                          "FileName"=paste(  dir.name, "/", classification, "/", classification, "-newick", "-", method, 
                                             ifelse( is.na(d.c.k$k[i]), "", paste( "=", sprintf("%.02f",d.c.k$k[i]), sep="") ), 
                                             ifelse( is.na(d.c.k$Distance[i]), "", paste( "+", d.c.k$Distance[i], sep="") ),
                                             ".csv", sep="" ),
                          stringsAsFactors=FALSE ));
      } );
      return (as.data.frame(do.call(rbind,ret.val)));
    } );
    all.methods <- as.data.frame(do.call(rbind, all.methods));
    
    # Read all the family trees from the topology file:
    cat( paste( "Caching the trees for classification '", classification, "'...\n", sep="") );
    roots <- languageclassification( classif.name=classification, csv.file=paste(dir.name,"/",classification,"/",classification,"-newick.csv",sep=""), csv.name.column="Family", csv.tree.column="Tree");
    if( is.null(roots) || length(roots) == 0 || !inherits(roots,"languageclassification") )
    {
      return (FALSE);
    }
    
    # Cache the trees for faster processing:
    roots.cached <- list();
    for( i in 1:nrow(all.methods) )
    {
      roots.cached[[i]] <- languageclassification( classif.name=all.methods$FileName[i], csv.file=all.methods$FileName[i], csv.name.column="Family", csv.tree.column="Tree" );
      names(roots.cached)[i] <- all.methods$FileName[i];
    }
    
    # All pairs of different methods:
    all.pairs.of.methods <- expand.grid( 1:nrow(all.methods), 1:nrow(all.methods) ); 
    all.pairs.of.methods <- all.pairs.of.methods[ all.pairs.of.methods[,1] > all.pairs.of.methods[,2], ]; # avoid comparisons with itself and allow only unique pairings
    all.pairs.of.methods <- all.pairs.of.methods[,c(2,1)];
    
    # For each pair of methods and family tree, compute the distances:
    cat( "    ... and computing the distances between the trees...\n" );
    ret.val <- mclapply( 1:nrow(all.pairs.of.methods), function(i)
    {
      #cat("i=",i,"\n");
      # Load the trees for this pair of methods:
      roots1 <- roots.cached[[ all.pairs.of.methods[i,1] ]];
      if( is.null(roots1) || length(roots1) == 0 || !inherits(roots1,"languageclassification") )
      {
        return (NULL);
      }
      roots2 <- roots.cached[[ all.pairs.of.methods[i,2] ]];
      if( is.null(roots2) || length(roots2) == 0 || !inherits(roots2,"languageclassification") )
      {
        return (NULL);
      }
      
      # For each family tree, compute the distances:
      ret.val <- lapply( 1:length(roots), function(j)
      {
        #cat( "j=",j,"\n");
        # Get the corresponding trees:
        tree.name <- get.name(roots[[j]]);
        tree1 <- find.language.family( roots1, tree.name );
        tree2 <- find.language.family( roots2, tree.name );
        # Compute the distances:
        d <- compare.trees(tree1, tree2);
        # Build the retun value:
        return (data.frame("Family"=tree.name,
                           "Classification"=classification,
                           "Method1"=all.methods$Method[all.pairs.of.methods[i,1]],
                           "Constant1"=all.methods$Constant[all.pairs.of.methods[i,1]],
                           "Distance1"=all.methods$Distance[all.pairs.of.methods[i,1]],
                           "Method2"=all.methods$Method[all.pairs.of.methods[i,2]],
                           "Constant2"=all.methods$Constant[all.pairs.of.methods[i,2]],
                           "Distance2"=all.methods$Distance[all.pairs.of.methods[i,2]],
                           "d.PH85"=d["PH85"], "d.score"=d["score"]));
      } );
      return (as.data.frame(do.call(rbind, ret.val)));
    }, mc.cores=parallel.mc.cores 
    );
    ret.val <- as.data.frame(do.call(rbind, ret.val));
    return( ret.val);  
  }
  compare.brlength.trees <- cmpfun(compare.brlength.trees);
  
  # Compare the trees for each of the classifications:
  for( classification in CLASSIFICATIONS )
  {
    # Apply only the relevant mg2015 distances:
    dists.codes <- DISTS.CODES[ -grep("mg2015", DISTS.CODES$Distance, fixed=TRUE), ]; dists.codes <- rbind(dists.codes, DISTS.CODES[ DISTS.CODES$Distance == paste0("mg2015(",classification,")"), ]);
    # Compare the trees for this classification:
    tmp <- compare.brlength.trees( dir.name="../output", classification=classification, 
                                   methods=METHODS, constant=CONSTANT, dists.codes=dists.codes,
                                   parallel.mc.cores=CPU_CORES, quotes=quotes );
    write.table( tmp, tree.comparisons, quote=FALSE, sep="\t", row.names=FALSE, col.names=(classification=="wals"), append=!(classification=="wals") ); # first classification only: overwrite the file
  }  
}






