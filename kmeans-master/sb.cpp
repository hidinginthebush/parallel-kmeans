#include <iostream>
#include  <omp.h>
#include<vector>
class qw{
    private:
        int a;
    public:
        qw(int n){
            a=n;
        }
        int get(){
            return a;
        }

};
int main()
{
	std::vector<qw> num;
    for(int i=0;i<10;i++){
        qw term(i);
        num.push_back(term);
    }
    num[0]=num[10];
    std::cout<<num[0].get();
    

}


