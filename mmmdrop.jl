#!/usr/bin/julia

### Specificaction
# mdrop submit-jobs
# mdrop attach-sessions
# mdrop push
# mdrop pull
# mdrop get-missing [folder]
# mdrop view-folder [folder]

include(homedir()*"/secrets")

CodeDir = "/home/janiserdmanis/Dropbox/Cebers/CodeBase/"
SESSION = "mmmdrop"

type PC
    sshlogin
    pasword
    workingdir
    maxcpuload
    actualload
    function PC(sshlogin,pasword,workingdir,maxcpuload)
        s = readall(`sshpass -p $pasword ssh -t $sshlogin "cat /proc/loadavg"`)
        actualload = parse(Float64,split(s," ")[3])
        println("pc $sshlogin has a load $actualload")
        new(sshlogin,pasword,workingdir,maxcpuload,actualload)
    end
end

pc1 = PC("janiserdmanis@5.179.6.146",PASS_PC1,"~/Dokumenti/calculation8/",6)
pc2 = PC("janiserdmanis@5.179.6.144",PASS_PC2,"~/Documents/calculation8/",4)

### Some convinient functions

function runpc(pc,command)
    readall(`sshpass -p $(pc.pasword) ssh -t $(pc.sshlogin) "$command"`)
end

function runtmux(pc,sesname,command)

    #run(`sshpass -p $(pc.pasword) ssh -t $(pc.sshlogin) "cd $(pc.workingdir)"`)
    
    try
        run(`sshpass -p $(pc.pasword) ssh -t $(pc.sshlogin) "tmux kill-session -t $sesname"`)
    catch
        info("session did not exist. Continuing with creating it")
    end
    run(`sshpass -p $(pc.pasword) ssh -t $(pc.sshlogin) "cd $(pc.workingdir) && tmux new-session -d -s $sesname && tmux send-keys -t $sesname:0 \"$command\" C-m "`)
end

function attachtotmux(pc,remoteses,localses,i)
    run(`tmux new-window -t $localses:$i`) # -n name
    run(`tmux send-keys -t $localses:$i "sshpass -p $(pc.pasword) ssh -t $(pc.sshlogin) \" tmux attach-session -t $remoteses \"" C-m`)
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

#replace(sesname,".","-")

### Setting up the commands to execute

function distribute(sesnames,commands;way="mirror")
    if way=="mirror"
        #for (i,comi) in enumerate(commands)
        for (sesn,comi) in zip(sesnames,commands)
            runtmux(pc1,sesn,comi)
            runtmux(pc2,sesn,comi)
        end
    elseif way=="cluster"
        for (sesn,comi) in zip(sesnames,commands)
            if pc1.actualload<pc1.maxcpuload
                info("Simulation $sesn submited on pc1")
                runtmux(pc1,sesn,comi)
                pc1.actualload += 1
            elseif pc2.actualload<pc2.maxcpuload
                info("Simulation $sesn submited on pc2")
                runtmux(pc2,sesn,comi)
                pc2.actualload += 1
            else
                info("Simulation $sesn can't be submited now")
            end
        end
    end
end

function collecttmux()
    try 
        run(`tmux kill-session -t $SESSION`)
    catch
        info("Local session was not found. Creating one.")
    end
    run(`tmux new-session -d -s $SESSION`)

    i = 1
    ses1 = getsessions(pc1)
    for ses in ses1
        attachtotmux(pc1,ses,SESSION,i)
        i+=1
    end

    ses2 = getsessions(pc2)
    for ses in ses2
        attachtotmux(pc2,ses,SESSION,i)
        i+=1
    end

    #run(`tmux attach-session -t $SESSION`)
end
    

#test()
#collecttmux()

### passing a folder with
function uniquecommands(loc)

    ses1 = getsessions(pc1)
    ses2 = getsessions(pc2)

    sess = []
    for ss in [ses1; ses2]

        if !(ss in sess)
            push!(sess,ss)
        else
            info("Dublicate session $ss")
        end
    end
    ### executes simulations from folder run if they are not present in the list sess
    ### executes the file manager.jl where as argument is given a specific simulation
    ### here one also can force to redo

    loc = normpath(loc)
    rr = search(loc,"calcpar/")
    relpath = loc[rr[1]:end]

    sesnames = []    
    commands = []
    for item in readdir(loc)
        if isfile(loc*item) && !(replace(item,".","-") in sess)
            comi = "julia manager.jl $relpath/$item"
            push!(commands,comi)
            push!(sesnames,replace(item,".","-"))
            #runtmux(pc1,item,comi)
            #runtmux(pc2,sesn,comi)
        end
    end
    return sesnames,commands
end

### truncating to a relative path
### I need to have a folder in the codebase because it stores a good exmaples and I also can avoid other means. 
### every configuration is stored to calcpar

# basedi = homedir()*"/Dropbox/Cebers/CodeBase/" ### I will get this from @__FILE__


### for debugging
# for i in 1:length(uniquecom)
#     uniquecom[i] = "echo "*uniquecom[i]
# end


# loc = normpath(loc)
# rr = search(loc,"calcpar/")
# relpath = loc[rr[1]:end]
# println(relpath)

#push()
#println(uniquecommands(loc))
#pull()

# distribute(sess,uniquecom,way="cluster")
# collecttmux()
#runtmux(pc1,"Bm35omega0.175.jl","ls")


function pushpc(pc)
    run(`sshpass -p $(pc.pasword) rsync -a --exclude="*/.*" --info=progress2  $CodeDir/ $(pc.sshlogin):$(pc.workingdir)`)
end

# http://serverfault.com/questions/219013/showing-total-progress-in-rsync-is-it-possible
function pullpc(pc)
    run(`sshpass -p $(pc.pasword) rsync -a --no-i-r --info=progress2 $(pc.sshlogin):~/SimulationData/ $(homedir())/SimulationData`)
end

if length(ARGS)>0
    if ARGS[1]=="attach-sessions"
        collecttmux()
    elseif ARGS[1]=="submit-jobs"
        loc = homedir()*"/Dropbox/Cebers/CodeBase/calcpar/run/"
        sess,uniquecom = uniquecommands(loc)
        distribute(sess,uniquecom,way="cluster")
    elseif ARGS[1]=="push"
        pushpc(pc1)
        pushpc(pc2)
    elseif ARGS[1]=="pull"
        pullpc(pc1)
        pullpc(pc2)
    elseif ARGS[1]=="test"
        sesnames = ["1","2","3"]
        commands = ["ls","pwd","ls"]
        distribute(sesnames,commands)
    elseif ARGS[1]=="clean-jobs" ### I will think about better way in the future
        #run(`echo rm -rf ~/Dropbox/Cebers/CodeBase/calcpar/run/*`)
        runpc(pc1,`rm -rf $(pc1.workingdir)/calcpar/run/*`)
        runpc(pc2,`rm -rf $(pc2.workingdir)/calcpar/run/*`)
    else
        error("No such command")
    end
end




