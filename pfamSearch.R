 pfamSearch <- function(
   
   fastas, 
   pfamA, 
   cut = 'ga', 
   cpus = 1L
   
   ){
  
  #Eval - err
  if (missing(fastas)){
    stop('No fasta files provided.')
  }
  
  if(length(fastas)<2){
    stop('At least 2 fasta files must be provided.')
  }
  
  if (missing(pfamA)){
    stop('Pfam-A.hmm file path is not provided.')
  }
  
  if (Sys.which("hmmsearch")==""){
    stop("\n\tHMMER (v.3) is not installed. (Couldn't find 'hmmsearch' in $PATH)
         \nPlease install it before re-running pfamSearch.\n\n")
  }

   
   #Load hmmpress wrapper function.
   hmmPress <- function (model) 
   {
     hmmpress <- paste("hmmpress -f", model)
     system(hmmpress, ignore.stdout = TRUE)
     o <- paste0(model, c("", ".h3f", ".h3i", ".h3m", ".h3p"))
     o
   }   

  #Load hmmsearch wrapper function.
  hmmSearch <- function (fasta, 
                         hmm, 
                         oty = "domtblout",
                         cut = "ga", 
                         n_threads = 1L) 
  {
    oty <- match.arg(oty, c("tblout", "domtblout"))
    cut <- match.arg(cut, c("ga", "tc"))
    blout <- paste0(sub('[.]faa$','',fasta), '_vs_pfamA')
    hmmse <- paste0("hmmsearch -o /dev/null --noali",
                    paste0(" --", oty, " "), 
                    blout, 
                    paste0(" --cut_", cut), 
                    paste0(" --cpu ",n_threads),
                    " ", 
                    hmm, 
                    " ", 
                    fasta)
    system(hmmse)
    return(blout)
  }
  
 
  
  #Press Pfam if not yet.
  idx <- paste0(pfamA, c('.h3f', '.h3i', '.h3m', '.h3p'))
  if (any(!file.exists(idx))){
    cat('Pfam-A.hmm is not indexed. Pressing Pfam-A.hmm.. ')
    hmmPress(model = pfamA)
    cat('DONE!\n')
  }
  
  #Hmmsearch
  cat('Searching.. ')
  domtblout <- parallel::mclapply(fastas, function(x){
    

    hmmres <- hmmSearch(fasta = x, 
                        hmm = pfamA, 
                        cut = cut, 
                        oty = 'domtblout', 
                        n_threads = 0)
    
    return(hmmres)
    
  }, mc.preschedule = FALSE, mc.cores = cpus)
 
  cat('DONE!\n')
  
  return(domtblout)
}


