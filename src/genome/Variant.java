package genome;

import java.util.HashMap;
import java.util.Random;

public class Variant {
	HashMap<String, String> _atts;
	private String _id, _ref;
	private int      _chr = -1, _pos = -1, _del = -1;
	private String[] _alts;
	private int      _maternal = 0, _paternal = 0;
	private boolean  _isPhased = false; // Phasing
	private static final Random _rand = new Random();

	
	public Variant(int chr,int pos,String ref,String[] alts,String phase){
		init(chr,pos,ref,alts,phase,new HashMap<String,String>(),"");
	}
	
	public Variant(int chr,int pos,String ref,String[] alts,String phase,HashMap<String, String> atts){
		init(chr,pos,ref,alts,phase,atts,"");
	}

	public Variant(int chr,int pos,String ref,String[] alts,String phase,HashMap<String, String> atts,String id){
		init(chr,pos,ref,alts,phase,atts,id);
	}

	public void init(int chr,int pos,String ref,String[] alts,String phase,HashMap<String, String> atts,String id){
		_atts = atts;
		_id   = id;
		_chr  = chr;
		_pos  = pos;
		_del  = ref.length();
		_ref  = ref;
		_alts = alts;

		phase = phase.trim();
		boolean strangePhase = false;
		if (phase.length() == 1) {
			if (Character.isDigit(phase.charAt(0))) {
				int val = Integer.parseInt(phase);
				if (chr == 22) {
					_maternal = val;
					_isPhased = true;
				} else if (chr == 23) {
					_paternal = val;
					_isPhased = true;
				} else strangePhase = true;
			} else strangePhase = true;
		} else if (phase.length() == 3) {
			char c1      = phase.charAt(0);
			char c2      = phase.charAt(2);
			char phasing = phase.charAt(1);

			if (Character.isDigit(c1) && Character.isDigit(c2)) {
				_paternal = Character.digit(c1,10);
				_maternal = Character.digit(c2,10);
				if (phasing == '|') _isPhased = true;
			} else strangePhase = true;
		} else strangePhase = true;

		if (strangePhase)
			System.err.println("Unreconized phasing '" + phase + "'.");
	}

	public String getID()		  { return _id;  }
	public int    getChromosome() { return _chr; }
	public int    getPosition()   { return _pos; }
	public int    getDeletion()   { return _del; }
	public String getReference()  { return _ref; }
	public int    getMaternal()   { return _maternal; }
	public int    getPaternal()   { return _paternal; }
	//	public String[] getAlternative(){ return _alts; } 
	public String getAlternative(){ 
		String ret=_alts[0]; 
		for(int i=1;i<_alts.length;i++)
			ret += _alts[i]+",";
		return ret;
	}

	public String getAttribute(String key){
		if (_atts.containsKey(key))
			return _atts.get(key);
		else
			return null;
	}
	
	
	
	public String insertion(int ind){
		if (ind <= 0 || ind > _alts.length) return "";
		return _alts[ind - 1];
	}
	public int variantBases() {
		int ret = _del;
		for (int i = 0;i < _alts.length;i++)
			if (_del != _alts[i].length()) ret += _alts[i].length();
		return ret;
	}
	public boolean isPhased() { return _isPhased; }
	public void randomizeHaplotype() {
		if (_rand.nextDouble() > 0.5) return;
		int tmp   = _paternal;
		_paternal = _maternal;
		_maternal = tmp;
		return;
	}
}
