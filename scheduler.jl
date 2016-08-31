### PCS are configured from the client side

# ssh-keygen
# ssh janiserdmanis@5.179.6.144; exit
# ssh-copy-id janiserdmanis@5.179.6.144
# And now it should be possibel to login without password

# In case of slow ssh connection times I added `UseDNS no` to /etc/ssh/sshd_config

type PC
    sshlogin
    workingdir
    maxcpuload
    actualload
    function PC(sshlogin,workingdir,maxcpuload)

        actualload = nothing

        s = @spawn readall(`ssh $sshlogin "cat /proc/loadavg"`)
        # wait(Timer(10))
        # if isready(s)
        s = fetch(s)
        actualload = parse(Float64,split(s," ")[1])
        # else
        #     info("Connection did not work")
        # end
        #println("pc $sshlogin has a load $actualload")
        new(sshlogin,workingdir,maxcpuload,actualload)
        ### return true ### or false depending if connection established succesfully.
    end
end

### have a dictionary of computers

pcs = Dict()

#pcs["UX305"] = PC("janiserdmanis@192.168.137.27",PASS_PC1,"~/Documents/calculation9/",4)
#pcs["JanisPC"] = PC("janiserdmanis@192.168.137.27","~/Documents/MDrop/",6)
pcs["JanisPC"] = PC("janiserdmanis@5.179.6.146","~/Dokumenti/MDrop/",6)
pcs["CimursPC"] = PC("janiserdmanis@5.179.6.144","~/Documents/MDrop/",3)

### See if command is already running

### Printing out all active simulations (for names using the previous convention), on which pc running

function runpc(pc,command)
    readall(`ssh $(pc.sshlogin) "$command"`)
end

function getsessions(pc)
    sesnames = Any[]

    s1 = runpc(pc,"tmux list-sessions 2>>/dev/null || true")
    ll = split(s1,"\n")
    ll = ll[1:length(ll)-1]
    for lli in ll
        se = search(lli,": ")[1]
        push!(sesnames,lli[1:se-1])
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

    run(`ssh -t $(pc.sshlogin) "cd $(pc.workingdir) && tmux new-session -d -s $sesname && tmux send-keys -t $sesname:0 \"$command\" C-m "`)
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

function parsesesname(command)
    command = replace(command," ","")
    command = replace(command,".","DOT")
end

function runtmux(command)
    sessions = getallsessions()

    sesnames = [sesi for (sesi,key) in sessions]
    if parsesesname(command) in sesnames
        error("Already running")
    else
        runtmux(parsesesname(command),command)
    end
end

### Printing out output of simulation -> connect to tmux
    
function pushpc(pc)
    run(`rsync -a --exclude="*/.*" --exclude=".git/*" --exclude="*/METADATA/*" --exclude="*/.git/*" --info=progress2 $BaseDir/ $(pc.sshlogin):$(pc.workingdir)`)
end

function pullpc(pc)
    run(`rsync -a --no-i-r --info=progress2 $(pc.sshlogin):~/SimulationData/ $(homedir())/SimulationData`)
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
        run(`ssh -t $(pc.sshlogin) "tmux attach-session -t $ses"`)    

    elseif ARGS[1]=="push" ### Additional argument for pushing last steps to all pcs
        for pc in values(pcs)
            pushpc(pc)
        end
    elseif ARGS[1]=="pull"
        for pc in values(pcs)
            pullpc(pc)
        end
    else
        error("Command $(ARGS[1]) not found.")
    end
    
end
   
