/****   
 *   -- treeSync --
 *   tool to synchronize two tag & probe trees so that they contain the same entries
 *   (e.g. for validation of re-reconstruction with different relases or conditions)
 *
 *   Takes as input 2 trees, produces as output 2 trees, which are copies of the
 *   input ones but preserving only the shared (tag+probe) pairs.
 *
 *   A pair P from the first is defined as shared if:
 *     1- the (run,lumi,event) are the same
 *     2- tag and probe (pt,eta,phi) are within some tolerance from each other
 *     3- this pair is the best match for all pairs from the second file
 *        that satisfy the first two criteria, and vice-versa.
 *   That is, if a pair P from the first file matches to two pairs P1, P2 in the second,
 *   only the best match is retained.
 *   (note: I didn't really prove that the algorithm converges for all possible inputs)
 *
 *  If you need to modify your defintion of pair, or the tolerances, modify the 'Item' class.
*/

///// Compile with
// gcc -O3 -std=c++0x -lstdc++ $(root-config --cflags --ldflags --libs) treeSync.cxx -o treeSync.exe
///// Run with
// ./treeSync.exe <dir> <file1> <file2>
//    where 
//        <dir> = directory inside of the files in which the 'probe_tree' tree is located
//        <file1>, <file2> = input files (must be local files, in a writable directory)
//    the output files will be named as <file1>, <file2> with '.root' -> '.shared.root'      
/////

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cassert>

/// utility: deltaPhi (unsigned)
inline float dphi(float phi1, float phi2) {
    float ret = std::abs(phi1-phi2);
    if (ret > M_PI) ret = 2*M_PI-ret;
    return ret;
}

/// Item for a tag+probe pair
struct Item { 
    static constexpr float tolerance = 0.1f;

    Item() : event(0) {}
    Item(unsigned int ev, float tagpt, float tageta, float tagphi, float probept, float probeeta, float probephi) :
        event(ev),
        tag_pt(tagpt), tag_eta(tageta), tag_phi(tagphi), probe_pt(probept), probe_eta(probeeta), probe_phi(probephi),
        count(1) {}

    unsigned int event; 
    float tag_pt, tag_eta, tag_phi, probe_pt, probe_eta, probe_phi; 
    mutable unsigned int count; 

    bool operator<(const Item &other) const {
        return event < other.event;
    }

    bool operator==(const Item &other) const {
        if (event != other.event) return false;
        if (std::abs(tag_eta   - other.tag_eta)   > tolerance) return false;
        if (std::abs(probe_eta - other.probe_eta) > tolerance) return false;
        if (dphi(tag_phi,other.tag_phi) > tolerance) return false;
        if (dphi(probe_phi,other.probe_phi) > tolerance) return false;
        if (std::abs(tag_pt - other.tag_pt)/(tag_pt + other.tag_pt) > tolerance) return false;
        if (std::abs(probe_pt - other.probe_pt)/(probe_pt + other.probe_pt) > tolerance) return false;
        return true;
    }
   
    float distance( const Item &other) const {
        if (event != other.event) return 999.0f;
        if (std::abs(tag_eta   - other.tag_eta)   > tolerance) return 99.0f;
        if (std::abs(probe_eta - other.probe_eta) > tolerance) return 99.0f;
        return  std::abs(tag_eta   - other.tag_eta) +
                std::abs(tag_eta   - other.tag_eta)   +
                std::abs(probe_eta - other.probe_eta) +
                dphi(tag_phi,other.tag_phi) +
                dphi(probe_phi,other.probe_phi) +
                std::abs(tag_pt - other.tag_pt)/(tag_pt + other.tag_pt) +
                std::abs(probe_pt - other.probe_pt)/(probe_pt + other.probe_pt);
    }

    void print(std::ostream& stream) const {
        using std::setw;
        using std::setprecision;
        using std::fixed;
        stream << " event " << setw(9) << event << " tag ( " << 
            setw(6) << fixed << setprecision(2) << tag_pt << ", " <<
            setw(5) << fixed << setprecision(2) << tag_eta << ", "  <<
            setw(5) << fixed << setprecision(2) << tag_phi << "), probe "  <<
            setw(6) << fixed << setprecision(2) << probe_pt << ", "  <<
            setw(5) << fixed << setprecision(2) << probe_eta << ", "  <<
            setw(5) << fixed << setprecision(2) << probe_phi << "), count  " << count;
    }
};
std::ostream& operator<<(std::ostream& stream, const Item& item) {
    item.print(stream);  
    return stream;
}

/// Class defining the interface for keeping count of tags and probes.
/// The actual implementation is in SlateImpl
template<class SlateImpl>
class Slate : public SlateImpl {
    public:
        typedef typename SlateImpl::value_type::value_type::iter_type  iter_type;
        typedef typename SlateImpl::value_type::value_type::range_type range_type;

        void add(unsigned int run, unsigned int ls, const Item item) {
            SlateImpl::get(run).get(ls).push_back(item);
        }

        const Item * get(unsigned int run, unsigned int ls, const Item &item) const { 
            if (!SlateImpl::has(run)) return 0;

            typename SlateImpl::value_type const & lumis = SlateImpl::get(run);
            if (!lumis.has(ls)) return 0;

            typename SlateImpl::value_type::value_type const & items = lumis.get(ls);
            return items.get(item);   
        }

        range_type equal_range(unsigned int run, unsigned int ls, const Item &item) const { 
            if (!SlateImpl::has(run)) return range_type();

            typename SlateImpl::value_type const & lumis = SlateImpl::get(run);
            if (!lumis.has(ls)) return range_type();

            typename SlateImpl::value_type::value_type const & items = lumis.get(ls);
            return items.equal_range(item);
        }
};

/// One implementation of a map from run to LumiList (which can be anything)
/// based on std::vector
template<class LumiList>
class VectorRunList {
    public:
        typedef LumiList value_type;
        VectorRunList() : offset_(0) {}
        VectorRunList(unsigned int min, unsigned int max) : 
            offset_(min), runs_(max-min+1) {}

        void init(unsigned int min, unsigned int max) { 
            offset_ = min; 
            runs_.clear(); 
            runs_.resize(max-min+1); 
        }

        LumiList & get(unsigned int run) {
            assert(has(run) && "Doesn't have the run!");
            return runs_[run-offset_]; 
        }
        const LumiList & get(unsigned int run) const { 
            assert(has(run) && "Doesn't have the run!");
            return runs_[run-offset_]; 
        }
        bool has(unsigned int run) const { 
            return  offset_ <= run && run < offset_+runs_.size(); 
        }

        void index() {
            typedef typename std::vector<LumiList>::iterator IT;
            for (IT it = runs_.begin(), ed = runs_.end(); it != ed; ++it)  {
                it->index();
            }
        }

    private:
        unsigned int offset_;
        std::vector<LumiList> runs_;
};

/// One implementation of a map from ls number to ItemList, based on std::vector
template<class ItemList>
class VectorLumiList {
    public:
        typedef ItemList value_type;
        VectorLumiList() : size_(0) {}
        ItemList & get(unsigned int ls) { 
            if (ls >= size_) { 
                lumis_.resize(std::max(ls+1, 2*size_)); 
                size_ = lumis_.size(); 
            }
            assert(ls < lumis_.size() && "Bogus lumi!");
            return lumis_[ls]; 
        }
        const ItemList & get(unsigned int ls) const { 
            assert(ls < lumis_.size() && "Bogus lumi!");
            return lumis_[ls]; 
        }
        bool has(unsigned int ls) const { return ls < size_; }

        void index() {
            typedef typename std::vector<ItemList>::iterator IT;
            for (IT it = lumis_.begin(), ed = lumis_.end(); it != ed; ++it)  {
                it->index();
            }
        }

    private:
        unsigned int size_;
        std::vector<ItemList> lumis_;
};

/// One implementation of list of items, based on std::vector
class SortedVectorItemList {
    public:
        typedef std::vector<Item>::const_iterator iter_type;
        typedef std::pair<iter_type,iter_type> range_type;
        SortedVectorItemList() {}

        /// range of items that _can_ be equal to item
        range_type equal_range(const Item &item) const {
            return std::equal_range(items_.begin(), items_.end(), item);
        }

        /// note: item pointer is volatile can be destroyed by any non-const op.
        const Item* get(const Item &item) const {
            range_type range = equal_range(item);
            if (range.first == range.second) return 0;
            iter_type best = range.first; float dmin = item.distance(*best); 
            for (++range.first; range.first != range.second; ++range.first) {
                float d =  item.distance(*range.first);
                if (d < dmin) { dmin = d; best = range.first; }
            }
            return (*best == item ? &*best : 0);
        }

        void push_back(const Item &item) {
            items_.push_back(item); 
            items_.back().count = 1; // important!
        }

        void index() { std::sort(items_.begin(), items_.end()); }
    protected:
        std::vector<Item> items_;
};

//// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TStopwatch.h>

void branches(TTree *t, unsigned int &run, unsigned int &lumi, Item &item) {
    item.count = 0;
    t->SetBranchStatus("*",0);
    t->SetBranchStatus("run",   1);
    t->SetBranchStatus("lumi",  1);
    t->SetBranchStatus("event", 1);
    t->SetBranchStatus("tag_pt",  1);
    t->SetBranchStatus("tag_eta", 1);
    t->SetBranchStatus("tag_phi", 1);
    t->SetBranchStatus("pt",  1);
    t->SetBranchStatus("eta", 1);
    t->SetBranchStatus("phi", 1);
    t->SetBranchAddress("run",   &run);
    t->SetBranchAddress("lumi",  &lumi);
    t->SetBranchAddress("event", &item.event);
    t->SetBranchAddress("tag_pt",  &item.tag_pt);
    t->SetBranchAddress("tag_eta", &item.tag_eta);
    t->SetBranchAddress("tag_phi", &item.tag_phi);
    t->SetBranchAddress("pt",  &item.probe_pt);
    t->SetBranchAddress("eta", &item.probe_eta);
    t->SetBranchAddress("phi", &item.probe_phi);
}

template<class AnySlate>
void fill(AnySlate &slate, TTree *t, std::vector<char> &ok) {
    unsigned int run, lumi; Item item;
    branches(t, run, lumi, item);

    unsigned int n = t->GetEntries();
    unsigned int rmin = t->GetMinimum("run");
    unsigned int rmax = t->GetMaximum("run");
    if (rmin < 0) rmin = 0; 
    if (rmax < 300000) rmax = 300000;
    slate.init(rmin,rmax);

    for (unsigned int i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (!ok[i]) continue;
        slate.add(run,lumi,item);
    }

    slate.index();
}

template<class AnySlate>
void sync(const AnySlate &slate, TTree *t, std::vector<char> &ok, int minCount=1) {
    unsigned int run, lumi; Item item;
    branches(t, run, lumi, item);

    unsigned int n = t->GetEntries();

    std::vector<char>::iterator itok = ok.begin();
    for (unsigned int i = 0; i < n; ++i, ++itok) {
        if (!*itok) continue;
        t->GetEntry(i);
        //if (run != 191226 || lumi != 237) continue;
        const Item *it = slate.get(run,lumi,item);
        if (it == 0 || it->count < minCount) { 
            *itok = false; 
        } else { 
            *itok = true; 
            it->count++; 
        }
    }
}

#include <TString.h>
void writeOut(TTree *tree, const std::vector<char> &ok, const char *fname, const char *dir) {
    TString outname = fname; outname.ReplaceAll(".root",".shared.root");
   
    TFile *fOut = new TFile(outname, "RECREATE");
    fOut->mkdir(dir)->cd();

    tree->SetBranchStatus("*",1); 
    TTree *tOut = tree->CloneTree(0);

    std::cout << "Writing output to " << fOut->GetName() << std::endl;
    int step = tree->GetEntries()/10;
    TStopwatch timer; timer.Start();
    for (int i = 0, n = tree->GetEntries(); i < n; ++i) {
        tree->GetEntry(i);
        if (ok[i]) tOut->Fill();
    }

    tOut->AutoSave(); // according to root tutorial this is the right thing to do

    fOut->Close();
}

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cerr << "Usage: treeSync dirName file1 file2 " << std::endl;
        return 1;
    }

    TFile *f1 = TFile::Open(argv[2]);
    TFile *f2 = TFile::Open(argv[3]);
    TTree *t1 = (TTree*)  f1->GetDirectory(argv[1])->Get("fitter_tree");
    TTree *t2 = (TTree*)  f2->GetDirectory(argv[1])->Get("fitter_tree");

    typedef Slate<VectorRunList<VectorLumiList<SortedVectorItemList> > > MySlate;
    MySlate slate;

    std::vector<char> ok1(t1->GetEntries(),1);
    std::vector<char> ok2(t2->GetEntries(),1);

    TStopwatch timer; timer.Start();
    for (;;) {
        unsigned int n1start = std::count(ok1.begin(),ok1.end(),true);
        unsigned int n2start = std::count(ok2.begin(),ok2.end(),true);
        std::cout << "First tree is from " << f1->GetName() << ": " << n1start << " entries" << std::endl;
        std::cout << "First tree is from " << f2->GetName() << ": " << n2start << " entries" << std::endl;

        //std::cout << "Filling in slate from first tree " << std::endl;
        fill(slate, t1, ok1);
        //std::cout << "Done after " << timer.RealTime()/60.0 << " min\n" << std::endl; timer.Continue();
    
        //std::cout << "Testing second tree vs slate " << std::endl;
        sync(slate, t2, ok2, 1);
        //std::cout << "Done after " << timer.RealTime()/60.0 << " min\n" << std::endl; timer.Continue();

        int nshared2 = std::count(ok2.begin(),ok2.end(),true);
        std::cout << "Entries of 2 matched at least to one entry of 1: " << nshared2 << std::endl; 

        //std::cout << "Testing first tree vs slate to see the match backs" << std::endl;
        sync(slate, t1, ok1, 2);
        //std::cout << "Done after " << timer.RealTime()/60.0 << " min\n" << std::endl; 

        int nshared1 = std::count(ok1.begin(),ok1.end(),true);
        std::cout << "Entries of 1 matched back from 2: " << nshared1 << std::endl; 

        if (nshared1 == nshared2) break;

        std::cout << "Swaping trees before repeating." << std::endl;
        std::swap(t1,t2);
        std::swap(f1,f2);
        std::swap(ok1,ok2);
    };

    std::cout << "Converged after " << std::setprecision(1) << std::fixed << timer.RealTime()/60.0 << " min." << std::endl;
    writeOut(t1, ok1, argv[2], argv[1]);
    writeOut(t2, ok2, argv[3], argv[1]);

    return 0;
}
