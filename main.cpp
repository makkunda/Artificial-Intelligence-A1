#include <iostream>
#include <string>
#include <cstdlib>
#include <sys/time.h>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

double proc_time;
int k,p,t; // k -> papers per session p -> no. of parallel sessions t-> no. of time intervals
double c;
int n;
vector<vector<double> > d;
string inputfilename,outputfilename;
double cur_score;
double best_score;
vector<int> curpos;
double MARGIN_TIME = 2000.0;	// in millis
timeval programStartTime, currentTime;
int fivepc;
int global_time = 0;
inline double millis(timeval &t) {
	return t.tv_sec*1e3 + t.tv_usec*1e-3;
}

double Temp ;

double compute_score(vector<int> &v)
{
    double score = 0.0;
    int i,j;
    for(i=0;i<n;i++)
    {
        int cur_time;

        cur_time = i/(p*k);

        int next_time_start = (cur_time + 1)*(p*k);

        int cur_session = (i%(p*k))/k;

        int next_session_start = cur_time*p*k + (cur_session+1)*k;

        for (j=(i+1);j<next_time_start;j++)
        {
            if(j<next_session_start)
            {
                score = score + 1 - d[v[i]][v[j]];
            }
            else
                score = score + c*d[v[i]][v[j]];
        }
    }
    return score;
}

void write_output()
{
    // if(cur_score>best_score)
    if(true)
    {
        ofstream fout(outputfilename.c_str());
        int i1,j1,k1;
        int cur = 0;
        for(i1=0;i1<t;i1++)
        {
            for(j1=0;j1<p;j1++)
            {
                for(k1=0;k1<k;k1++)
                {
                    fout<<curpos[cur]<<" ";
                    cur++;
                }
                if(j1 <(p-1))
                    fout<<"| ";
            }
            fout<<"\n";
        }
        best_score = cur_score;
        double now_score = compute_score(curpos);
        cout<<"now score "<<now_score;
    }
}


double change_score_swap(int i,int j,vector<int> &v)
{
    int time_i,time_j,sess_i,sess_j;
    // cout<<"before times computed for swap "<<endl;
    time_i = i/(p*k);
    time_j = j/(p*k);
    // cout<<"times computed for swap "<<endl;
    sess_i = (i%(p*k))/k;
    sess_j = (j%(p*k))/k;
    // cout<<"sess computed for swap "<<endl;
    int i1,j1;
    double ci1=0.0,ci2=0.0,cj1=0.0,cj2=0.0;

    for(i1=(time_i*p*k);i1<((time_i+1)*p*k);i1++)
    {
        if(i1!=i)
        {
            int cur_sess = (i1%(p*k))/k;
            if(cur_sess==sess_i)
            {
                ci1 = ci1 + 1 - d[v[i]][v[i1]];
                cj2 = cj2 + 1 - d[v[j]][v[i1]];
            }
            else
            {
                ci1 = ci1 +  c*d[v[i]][v[i1]];
                cj2 = cj2 +  c*d[v[j]][v[i1]];
            }
        }        
    }

    for(i1=(time_j*p*k);i1<((time_j+1)*p*k);i1++)
    {
        if(i1!=i)
        {
            int cur_sess = (i1%(p*k))/k;
            if(cur_sess==sess_j)
            {
                ci2 = ci2 + 1 - d[v[i]][v[i1]];
                cj1 = cj1 + 1 - d[v[j]][v[i1]];
            }
            else
            {
                ci2 = ci2 +  c*d[v[i]][v[i1]];
                cj1 = cj1 +  c*d[v[j]][v[i1]];
            }
        }        
    }

    double delta = ci2 + cj2 - ci1 - cj1;

    return delta;
}


double restart_shuff(vector<int> &v,int shuffle_size)
{
 int i1,j1,k1;
 double tot_change = 0;
 for(i1=0;i1<shuffle_size;i1++)
 {
     int rand1,rand2;
     // compute rand
     rand1 = rand();
     rand2 = rand();
     j1 = rand1 % n;
    if (j1 > 0)
        k1 = rand2 % j1;
    else
        k1 = rand2 % n;

     int temp1 = v[j1];
     v[j1] = v[k1];
     v[k1] = temp1;

     tot_change = tot_change + change_score_swap(j1,k1,v);
 }
 return tot_change;
}


void step(double stop_time)
{
    // cout<<"stop time is "<<stop_time<<endl;
    timeval diffTime;
    bool cont = true;
    bool restart = false;
    bool change = false;
    while(cont)
    {
        // cout<<"inside while loop "<<global_time<<endl;
        global_time ++;
        if(restart)
        {
            global_time = global_time/2;
            if(change)
                write_output();
            cur_score = cur_score + restart_shuff(curpos,fivepc);
            restart = false;
            change = false;
        }
        else
        {
            int try_count = 100;
            restart = true;
            int i1;
            for(i1=0;i1<try_count;i1++)
            {
                int j1,k1;
                int rand1,rand2;
                // compute rand
                rand1 = rand();
                rand2 = rand();

                // cout<<"rand computed for time "<<i1<<endl;
                j1 = rand1 % n;
                if (j1 > 0)
                    k1 = rand2 % j1;
                else
                    k1 = rand2 % n;
                // cout<<"rand modded for time "<<i1<<endl;
                double delta = change_score_swap(j1,k1,curpos);
                // cout<<"delta computed"<<endl;
                bool do_swap ;

                if(delta > 0)
                    do_swap = true;
                
                else{
                    // cout<<"before simul anneal"<<endl;
                    double cur_temp = Temp/((double)global_time);
                    double prob = exp(delta/cur_temp);

                    double sample = ((double)rand())/((double)RAND_MAX);

                    if (prob<sample)
                        do_swap = true;
                    else
                        do_swap = false;
                    // cout<<"after simul anneal"<<endl;
                }

                if(do_swap)
                {
                    int temp1 = curpos[j1];
                    curpos[j1] = curpos[k1];
                    curpos[k1] = temp1;

                    cur_score = cur_score + delta;
                    restart = false;
                    change = true;
                    break;
                }
            }
        }
        gettimeofday(&currentTime, NULL);
		timersub(&currentTime, &programStartTime, &diffTime);
        if(millis(diffTime) >= stop_time - MARGIN_TIME)
            cont = false;
    }
    write_output();
    cout<<cur_score<<endl;
}



int main ( int argc, char** argv )
{
    // cout<<"entered main"<<endl;
    gettimeofday(&programStartTime, NULL);
    // cout<<"read time"<<endl;
	srand(time(NULL));
    // cout<<"rand initialized"<<endl;
    inputfilename =  argv[1];
    outputfilename = argv[2];
    // cout<<"input arguments setup"<<endl;
    ifstream fin(inputfilename.c_str());
    
    fin>>proc_time;
    proc_time = proc_time*60000;
    fin>>k>>p>>t;
    fin>>c;
    n = k*p*t;
    d.resize(n);
    int i,j;
    for(i=0;i<n;i++)
        d[i].resize(n);
    
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
            fin>>d[i][j];
    }
    // cout<<"read input"<<endl;
    fivepc = (int)(0.05*((double)n));
    fivepc = max(5,fivepc);
    
    // set temp for simulated annealing
    Temp = 400.0;

    curpos.resize(n);
    for(i=0;i<n;i++)
        curpos[i] = i;
    cur_score = compute_score(curpos);
    best_score = (double)INT16_MIN;
    // cout<<"before first call"<<endl;
    step((proc_time/2.0));
    step(proc_time);
    return 0;
}

