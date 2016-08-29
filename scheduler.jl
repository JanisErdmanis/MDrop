include(homedir()*"/secrets")

type PC
    sshlogin
    pasword
    workingdir
    maxcpuload
    actualload
    function PC(sshlogin,pasword,workingdir,maxcpuload)

        actualload = nothing

        s = @spawn readall(`sshpass -p $pasword ssh $sshlogin "cat /proc/loadavg"`)
        # wait(Timer(10))
        # if isready(s)
        s = fetch(s)
        actualload = parse(Float64,split(s," ")[3])
        # else
        #     info("Connection did not work")
        # end
        #println("pc $sshlogin has a load $actualload")
        new(sshlogin,pasword,workingdir,maxcpuload,actualload)
        ### return true ### or false depending if connection established succesfully.
    end
end


### have a dictionary of computers

pcs = Dict()

#pcs["UX305"] = PC("janiserdmanis@192.168.137.27",PASS_PC1,"~/Documents/calculation9/",4)
#pcs["JanisPC"] = PC("janiserdmanis@192.168.137.27",PASS_PC1,"~/Documents/calculation9/",6)
pcs["JanisPC"] = PC("janiserdmanis@5.179.6.146",PASS_PC1,"~/Dokumenti/calculation9/",6)
pcs["CimursPC"] = PC("janiserdmanis@5.179.6.144",PASS_PC2,"~/Documents/calculation9/",4)

### See if command is already running

### Printing out all active simulations (for names using the previous convention), on which pc running

function runpc(pc,command)
    readall(`sshpass -p $(pc.pasword) ssh $(pc.sshlogin) "$command"`)
end

function getsessions(pc)
    sesnames = Any[]
    try
        s1 = runpc(pc,"tmux list-sessions")
        ll = split(s1,"\n")
        ll = ll[1:length(ll)-1]
        for lli in ll
            se = search(lli,": ")[1]
            push!(sesnames,lli[1:se-1])
        end
    end
    return sesnames
end

function getallsessions()

    allsessions = Any[]
    
    for (key,pc) in pcs
        if pc.actualload==nothing
            continue
        end
        sessions = getsessions(pc)
        for j in 1:length(sessions)
            push!(allsessions,("$(sessions[j])","$key"))
        end
    end

    return allsessions
end
    

### Printing out load for pcs

### Executting command on remote pc in background

function runtmux(pc,sesname,command)

    run(`sshpass -p $(pc.pasword) ssh -t $(pc.sshlogin) "cd $(pc.workingdir) && tmux new-session -d -s $sesname && tmux send-keys -t $sesname:0 \"$command\" C-m "`)
end

function runtmux(sesname,command)

    # Function which checks if sesname is not already running

    submitted = false

    for (key,pc) in pcs
        if pc.actualload<pc.maxcpuload
            runtmux(pc,sesname,command)
            info("Simulation $sesname submited on $key")
            submitted = true
            break
        end
    end

    if submitted==false
        error("PCS are overloaded with tasks")
    end
end

function runtmux(command)
    sessions = getallsessions()

    sensnames = [sesi for (sesi,key) in sessions]
    if command in sesnames
        error("Already running")
    else
        runtmux(command,command)
    end
end

### Printing out output of simulation -> connect to tmux
    
function pushpc(pc)
    run(`sshpass -p $(pc.pasword) rsync -a --exclude="*/.*" --info=progress2  $BaseDir/ $(pc.sshlogin):$(pc.workingdir)`)
end

function pullpc(pc)
    run(`sshpass -p $(pc.pasword) rsync -a --no-i-r --info=progress2 $(pc.sshlogin):~/SimulationData/ $(homedir())/SimulationData`)
end

if length(ARGS)==0
    ### Prints status by default
    println("###### CPU load of available pcs #######")

    for (key,pc) in pcs
        println("$key $(pc.sshlogin) has a load $(pc.actualload)")
    end
    
    println("###### Active calculations #######")

    allsesnames = getallsessions()

    for (i,(sesi,key)) in enumerate(allsesnames)
        println("$i \t $sesi \t $key")
    end

else
    # !(replace(item,".","-") 
    
    if ARGS[1]=="execute"
        if length(ARGS)==1
            error("no argument")
        else
            runtmux(ARGS[2])
        end
        
    elseif ARGS[1]=="connect"

        sessions = getallsessions()
        if length(ARGS)>1
            choice = parse(Int,ARGS[2])
        else
            for (i,(sesi,key)) in enumerate(sessions)
                println("$i \t $sesi \t $key")
            end

            choice = parse(Int,readline(STDIN))
        end
        ses,key = sessions[choice]
        pc = pcs[key]
        run(`sshpass -p $(pc.pasword) ssh -t $(pc.sshlogin) "tmux attach-session -t $ses"`)    

    elseif ARGS[1]=="push" ### Additional argument for pushing last steps to all pcs
        for pc in values(pcs)
            pushpc(pc)
        end
    elseif ARGS[1]=="pull"
        for pc in values(pcs)
            pullpc(pc)
        end
    end
end
   
