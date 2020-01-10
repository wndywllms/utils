import subprocess
import os
import pwd
import time

scriptpath = '/home/{user}/para/scripts/'.format(user=pwd.getpwuid(os.getuid())[0])

def run_parallel(cmds, logs='none', Monitor=False, monitorname="run", nthreads=2, sleeptime=0.5):
    '''run a list of shell commands in run_parallel
optional - specify max number of processes (default 2)
           secify file to log output to for each command
    '''

    processes = set()
    max_processes = nthreads
    process_monitors = set()


    for icmd, cmd in enumerate(cmds):
        if logs not in ['auto', 'none']:
            cmdlog = logs[icmd]
        elif logs == 'auto':
            cmdlog = "log{icmd:d}".format(icmd=icmd)
        else:
            cmdlog = ''
        if logs != 'none':
            with open(cmdlog,'w') as f:
                print('exec {icmd:d}: {cmd} > {log}'.format(cmd=cmd, icmd=icmd, log=cmdlog))
                p = subprocess.Popen(cmd.split(), stdout=f, stderr=subprocess.STDOUT)
                processes.add(p)
                pid=p.pid
                time.sleep(sleeptime)
                
                if Monitor:
                    monitor_cmd = 'python {scriptpath}/monitorjob.py {pid} {monitorname}'.format(pid=pid, monitorname=monitorname, scriptpath=scriptpath)
                    process_monitors.add(subprocess.Popen(monitor_cmd.split()))
        else:
            print('exec {icmd:d}: {cmd} > {log}'.format(cmd=cmd, icmd=icmd, log=cmdlog))
            p = subprocess.Popen(cmd.split())
            time.sleep(sleeptime)
            processes.add(p)
            pid=p.pid
                
            if Monitor:
                monitor_cmd = 'python {scriptpath}/monitorjob.py {pid} {monitorname}'.format(pid=pid, monitorname=monitorname, scriptpath=scriptpath)
                process_monitors.add(subprocess.Popen(monitor_cmd.split()))
        if len(processes) >= max_processes:
            os.wait()
            #if logs != 'none':
                #for ip, p in enumerate(processes):
                    #cmdlog = logs[ip]
                    #with open(cmdlog,'w') as f:
                        #ss = p.communicate()[0]
                        #try:
                            #f.write(ss)
                        #except:
                            #print ss
                            #print "log writing failed"
            processes.difference_update([p for p in processes if p.poll() is not None])
    for p in processes:
        if p.poll() is None:
            p.wait()
    # if we have submitted < max_processes, or on the last set #
    #os.wait()
    #if logs is not None:
        #for ip,p in enumerate(processes):
            #cmdlog = logs[ip]
            #with open(cmdlog,'w') as f:
                #f.write(p.communicate()[0])
    #processes.difference_update([p for p in processes if p.poll() is not None])
            
    return


def run_os_command(cmd):
    '''check for failure
    '''
    logging.info(cmd)
    
    out = os.system(cmd)
    if out != 0:
        print("ERROR")
        sys.exit()
    return