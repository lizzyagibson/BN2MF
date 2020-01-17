# sampler loop
# This version records everything
#
# Harness :: SeedOverRight , InitState :: () -> NewState , transitionProposal :: (State-> StateChange)
#     , ApplytransitionUpdate :: (CurrentState, Proposal) -> NewState
#     , shouldItTerminate :: ( IterationCount,CurrentState,CurrentProposal -> Yes or NO  )
#  -> (finalState, Number of Steps, Initial rng state,Initial State Value, stateHistory, transitionProposalHistory)

harness <- function(seed=NA # seed is fine to leave as NA
    ,InitState=NA
    ,TransitionProposal=NA
    ,ApplyTransition=NA
    ,ShouldWeTerminate=NA ) {
suppressWarnings({

if(is.na(InitState)
    || is.na(TransitionProposal)
    || is.na(ApplyTransition)
    || is.na(ShouldWeTerminate)){
  errorCondition(message=
    "you must use:
      InitState
      TransitionProposal
      ApplyTransition
      ShouldWeTerminate
    in keyword argument format
    "
    )
}
initialRngSnapshot <- tryCatch({.Random.seed},error=function(e){rnorm(1); .Random.seed})

restoreRNG <- FALSE

if (!is.na(seed))
  { set.seed(seed)
    restoreRNG <- TRUE
  }else {seed = initialRngSnapshot}

}) # end warning suppression

CurrentStep <- 0
CurrentState <- InitState()
StateRecord <- c(CurrentState) # this adds each state record, so keeps all
CurrentProposal <- TransitionProposal(CurrentState) # this does the first update step
ProposalRecord <- c(CurrentProposal) # This puts the output of the last line as the first record in a new list

while(!ShouldWeTerminate(CurrentStep,CurrentState,CurrentProposal)){
  CurrentStep <- CurrentStep + 1
  CurrentState <- ApplyTransition(CurrentState,CurrentProposal) # Update parameters
  StateRecord <- c(StateRecord,CurrentState) # Add current state to state record list
  CurrentProposal <- TransitionProposal(CurrentState) # Gibbs always accepts transitions
  ProposalRecord <- c(ProposalRecord,CurrentProposal) # Add current proposal to record list
  # second to last record proposed is final record
}

# restore RNG if required
if(restoreRNG)
  {set.seed(initialRngSnapshot)}

list(FinalState=CurrentState,StepCount=CurrentStep,TheStateRecord=StateRecord,TheProposalRecord=ProposalRecord,InitialRNGState=seed)
}

## Trivial example

# Gauss1dimWalk <- function(){
#   harness(
#     ,InitState=function(){rnorm(1)}
#     ,TransitionProposal=function(x){rnorm(1) }# doesn't care about current state
#        # for graphs, you might depend on the # of nbrs of current vertex/state to determine how you sample over choosing which neighbor you propose moving to
#     ,ApplyTransition=function(state,proposal){state + proposal }
#     ,ShouldWeTerminate=function(step,state,proposal){ (step > 10 )|| sum(abs(proposal)) > 100 || sqrt(sum(state^2) > 400)
#     }
#     )
# }
# 
# Gauss1dimWalk()
