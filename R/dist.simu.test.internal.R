dist.simu.test.internal <-
function(X, S1, S2, S1args, S2args, index)
    {
    ICS <- ics2(X, S1 = S1, S2 = S2, S1args = S1args, S2args = S2args)
    ics.distances(ICS, index = index)
    }
