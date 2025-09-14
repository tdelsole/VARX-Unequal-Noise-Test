sync.periods = function(dates1,dates2) {
  ### FINDS THE DATES THAT ARE COMMON BETWEEN TWO DATE VECTORS
  ### FOR MONTHLY DATA, 'DAY' SHOULD BE A COMMON VALUE (E.G., 1 OR 15)
  ### INPUT:
  #     dates1[ntot1]: first set of dates.  
  #     dates2[ntot2]: second set of dates
  ### OUTPUT: LIST
  #     $dates.common[ndates.common]: dates that are common
  #     $ndates.common: length of dates.common
  #     $nst1: index identifying the beginning of common period in dates1
  #     $nnd1: index identifying the end of       common period in dates1  
  #     $nst2: index identifying the beginning of common period in dates2
  #     $nnd2: index identifying the end of       common period in dates2
  
  dates1 = as.Date(dates1); dates2 = as.Date(dates2)
  day(dates1) = 15; day(dates2) = 15
  
  dates.common = as.Date(intersect(dates1,dates2))
  # if (any( diff(dates.common) < 28 | diff(dates.common) > 31)) stop('dates.common not consecutive')
  if (diff(range(diff(dates.common))) > 3) stop('dates.common not consecutive')
  
  ndates.common = length(dates.common)
  nst1       = which(dates.common[             1] == dates1)
  nnd1       = which(dates.common[ndates.common]  == dates1)
  nst2       = which(dates.common[            1]  == dates2)
  nnd2       = which(dates.common[ndates.common]  == dates2)
  sync.check = dates1[nst1:nnd1] == dates2[nst2:nnd2]
  if (any(!sync.check)) stop('could not synchronize the data')
  
  list(dates.common = dates.common, ndates.common = ndates.common, 
       nst1 = nst1, nnd1 = nnd1, nst2 = nst2, nnd2 = nnd2)
  
}