// this function is to try to track the Cs clocks in gps network to eliminate them from network
// created by tdaykin 09/29/20

int cs_count;
cs_count = 0;
std::vector<std::string> tmpsvn;
std::vector<double> tmpsvn;
for(i=0;i<num_clocks;i++){
    if(clk[i]=='Cs'){
        tmpsvn[i] = svn[i];
        cs_count++;
        std::cout << "this is SVN of Cs clock: " << tmpsvn[i] << "\n";
    }
}



