include(homedir()*"/secrets")

type PC
    sshlogin
    pasword
    workingdir
    maxcpuload
    actualload
    function PC(sshlogin,pasword,workingdir,maxcpuload)
        s = readall(`sshpass -p $pasword ssh $sshlogin "cat /proc/loadavg"`)
        actualload = parse(Float64,split(s," ")[3])
        #println("pc $sshlogin has a load $actualload")
        new(sshlogin,pasword,workingdir,maxcpuload,actualload)
    end
end

pc0 = PC("janiserdmanis@192.168.137.27",PASS_PC1,"~/Documents/calculation9/",4)
# pc1 = PC("janiserdmanis@5.179.6.146",PASS_PC1,"~/Dokumenti/calculation9/",6)
# pc2 = PC("janiserdmanis@5.179.6.144",PASS_PC2,"~/Documents/calculation9/",4)

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

### Printing out load for pcs

### Executting command on remote pc in background

function runtmux(pc,sesname,command)

    #run(`sshpass -p $(pc.pasword) ssh -t $(pc.sshlogin) "cd $(pc.workingdir)"`)
    
    try
        run(`sshpass -p $(pc.pasword) ssh -t $(pc.sshlogin) "tmux kill-session -t $sesname"`)
    catch
        info("session did not exist. Continuing with creating it")
    end
    run(`sshpass -p $(pc.pasword) ssh -t $(pc.sshlogin) "cd $(pc.workingdir) && tmux new-session -d -s $sesname && tmux send-keys -t $sesname:0 \"$command\" C-m "`)
end

function runtmux(sesname,command)

    ### Function which checks if sesname is not already running
    # error("Command already runing")

    if pc0.actualload<pc0.maxcpuload
        info("Simulation $sesname submited on pc0")
        runtmux(pc0,sesname,command)
        pc0.actualload += 1
    # elseif pc1.actualload<pc1.maxcpuload
    #     info("Simulation $sesname submited on pc1")
    #     runtmux(pc1,sesname,command)
    #     pc1.actualload += 1
    # elseif pc2.actualload<pc2.maxcpuload
    #     info("Simulation $sesname submited on pc2")
    #     runtmux(pc2,sesname,command)
    #     pc2.actualload += 1
    else
        error("PCS are overloaded with tasks")
    end
end

function runtmux(command)
    sesnames = getsessions(pc0)

    if command in sesnames
        error("Already running")
    else
        runtmux(command,command)
    end
end

# function runtmux(command)

# end

### Printing out output of simulation -> connect to tmux

function connectsession()

    sessions = getsessions(pc0)

    for i in 1:length(sessions)
        println("$i \t $(sessions[i]) \t pc0")
    end

    choice = parse(Int,readline(STDIN))
    pc = pc0
    run(`sshpass -p $(pc.pasword) ssh -t $(pc.sshlogin) "tmux attach-session -t $(sessions[choice])"`)    
end
    
function pushpc(pc)
    run(`sshpass -p $(pc.pasword) rsync -a --exclude="*/.*" --info=progress2  $BaseDir/ $(pc.sshlogin):$(pc.workingdir)`)
end

function pullpc(pc)
    run(`sshpass -p $(pc.pasword) rsync -a --no-i-r --info=progress2 $(pc.sshlogin):~/SimulationData/ $(homedir())/SimulationData`)
end

if length(ARGS)==0
    ### Prints status by default
    println("###### CPU load of available pcs #######")

    pc = pc0
    println("pc $(pc.sshlogin) has a load $(pc.actualload)")
    
    println("###### Active calculations #######")
    sessions = getsessions(pc0)

    for i in 1:length(sessions)
        println("$i \t $(sessions[i]) \t pc0")
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

        sessions = getsessions(pc0)
        if length(ARGS)>1
            choice = parse(Int,ARGS[2])
        else
            for i in 1:length(sessions)
                println("$i \t $(sessions[i]) \t pc0")
            end

            choice = parse(Int,readline(STDIN))
        end
        pc = pc0
        run(`sshpass -p $(pc.pasword) ssh -t $(pc.sshlogin) "tmux attach-session -t $(sessions[choice])"`)    

    elseif ARGS[1]=="push" ### Additional argument for pushing last steps to all pcs
        pushpc(pc0)
    elseif ARGS[1]=="pull"
        pullpc(pc0)
    end
end
    
