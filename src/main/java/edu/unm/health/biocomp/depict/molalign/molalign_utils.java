package edu.unm.health.biocomp.depict.molalign;

import java.io.*;
import java.util.*;
import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.util.*; //MolHandler,MolAligner
import chemaxon.sss.search.*; //MolSearch, SearchException

import com.chemaxon.search.mcs.*; //MaxCommonSubstructure,McsSearchOptions,McsSearchResult

/////////////////////////////////////////////////////////////////////////////
/**	Static methods for molecular alignment.
	NOTE: MCES is gone (between 5.8.3 and 6.3.1). Now MaxCommonSubstructure
	@author Jeremy J Yang
*/
public class molalign_utils
{
  /////////////////////////////////////////////////////////////////////////////
  public static boolean AlignToMCS(MaxCommonSubstructure mcs,McsSearchResult mcs_result)
  {
    Molecule qmol = mcs.getQuery();
    Molecule tmol = mcs.getTarget();
    if (qmol.getDim()!=2) qmol.clean(2,null,null);
    if (tmol.getDim()!=2) tmol.clean(2,null,null);
    MolAligner molaligner = new MolAligner();
    molaligner.setPatternMolecule(qmol);
    molaligner.setTargetMolecule(tmol);
    int [] amap = mcs_result.getAtomMapping();
    if (amap==null) return false;
    molaligner.align(amap);
    return true;
  }
  /////////////////////////////////////////////////////////////////////////////
  public static int AlignToMCES(ArrayList<Molecule> mols)
  {
    McsSearchOptions mcsOpts = new McsSearchOptions.Builder()
	.connectedMode(true)
	.minFragmentSize(5)
	.build();
    MaxCommonSubstructure mcs = MaxCommonSubstructure.newInstance(mcsOpts);

    mcs.setQuery(mols.get(0));
    int n_aligned=0;
    for (int i=1;i<mols.size();++i)
    {
      mcs.setTarget(mols.get(i));
      McsSearchResult mcs_result = mcs.nextResult();
      if (AlignToMCS(mcs,mcs_result)) ++n_aligned;
    }
    return n_aligned;
  }
  /////////////////////////////////////////////////////////////////////////////
  public static int AlignToSmarts(ArrayList<Molecule> mols,String smarts)
  {
    MolHandler smartsReader = new MolHandler();
    smartsReader.setQueryMode(true);
    try { smartsReader.setMolecule(smarts); }
    catch (MolFormatException e) {
      System.err.println(e.getMessage());
      return 0;
    }
    Molecule qmol = smartsReader.getMolecule();
    MolSearch search = new MolSearch();
    search.setQuery(qmol);
    search.setTarget(mols.get(0));
    int [] hitatoms=null;
    try { hitatoms=search.findFirst(); }
    catch (SearchException e) {
      System.err.println(e.getMessage());
      return 0;
    }
    if (hitatoms==null) { return 0; } // no match, no aligning.
    qmol.clean(2,null,null);
    for (int i=0;i<qmol.getAtomCount();++i)
    {
      qmol.getAtom(i).setX(mols.get(0).getAtom(hitatoms[i]).getX());
      qmol.getAtom(i).setY(mols.get(0).getAtom(hitatoms[i]).getY());
      qmol.getAtom(i).setZ(mols.get(0).getAtom(hitatoms[i]).getZ());
    }
    MolAligner molaligner = new MolAligner();
    molaligner.setPatternMolecule(qmol);
    int n_aligned=0;
    for (Molecule mol: mols)
    {
      search.setTarget(mol);
      try { hitatoms=search.findFirst(); }
      catch (SearchException e) { continue; }
      if (hitatoms==null) { continue; } // no match, no aligning.
      molaligner.setTargetMolecule(mol);
      molaligner.align(hitatoms); // alternative (1)

      //molaligner.setTargetMolecule(mol); // alternative (2)
      //MolAligner.AlignmentResult alignment = molaligner.calculate(hitatoms); // alternative (2)
      //molaligner.rotate(alignment); // alternative (2)

      ++n_aligned;
    }
    return n_aligned;
  }
}
