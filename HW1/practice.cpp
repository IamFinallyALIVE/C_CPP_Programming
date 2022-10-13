#include<iostream>
#include<map>
using namespace std;


int main(){


    map<int,int>quiz;

    quiz.insert(pair<int,int>(1,40));
    quiz.insert(pair<int,int>(2,0));
    quiz.insert(pair<int,int>(4,30));

    auto it = quiz.find(9);
    cout<<it->second<<endl;
    map<int,int>::iterator itr;

    for(itr = quiz.begin();itr!=quiz.end();itr++){

        cout<<itr->first<<" "<<itr->second<<endl;
    }
    return 0;
}