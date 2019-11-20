package edu.unm.health.biocomp.depict;

import java.io.*;
import java.lang.Math;
import java.net.*; //URLEncoder,InetAddress
import java.text.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.*;
import java.awt.Color;
import javax.servlet.*;
import javax.servlet.http.*;

import com.oreilly.servlet.MultipartRequest;
import com.oreilly.servlet.multipart.DefaultFileRenamePolicy;
import com.oreilly.servlet.*; //Base64Encoder,Base64Decoder

import chemaxon.formats.*;
import chemaxon.util.MolHandler;
import chemaxon.struc.*;
import chemaxon.struc.prop.MMoleculeProp;
import chemaxon.sss.search.*; //MolSearch,SearchException
import chemaxon.license.*; //LicenseException,LicenseManager
import chemaxon.marvin.io.MolExportException;

import edu.unm.health.biocomp.util.*;
import edu.unm.health.biocomp.util.http.*;
import edu.unm.health.biocomp.depict.molalign.*;

/**	Depict molecules as PNG or JPEG images using mol2img servlet.

	Capabilities:<UL>
	<LI>  Use input 2D coords via mol2img (various 2D formats)
	<LI>  Handle depiction of smarts (mol2img can do it).
	<LI>  Transparent PNG.
	<LI>  MRV abbreviations (e.g. OEt, CF3) (Note: -H side-effect disables this).
	<LI>  Align to smarts match.
	</UL>
	@author Jeremy J Yang
*/
public class depict_servlet extends HttpServlet
{
  private static String SERVLETNAME=null;
  private static String CONTEXTPATH=null;
  private static String LOGDIR=null;	// configured in web.xml
  private static String APPNAME=null;	// configured in web.xml
  private static String UPLOADDIR=null;	// configured in web.xml
  private static int N_MAX=100; // configured in web.xml
  private static ServletContext CONTEXT=null;
  private static ResourceBundle rb=null;
  private static PrintWriter out=null;
  private static ArrayList<String> outputs=null;
  private static ArrayList<String> errors=null;
  private static HttpParams params=null;
  private static ArrayList<Molecule> mols=null;
  private static LinkedHashMap<String,Integer> sizes_h=null;
  private static LinkedHashMap<String,Integer> sizes_w=null;
  //private static Integer SERVERPORT=null;
  private static String SERVERNAME=null;
  private static String REMOTEHOST=null;
  private static String DATESTR=null;
  private static File logfile=null;
  private static String color1="#EEEEEE";
  private static MolSearch molsearch=null;
  private static ArrayList<Color> atomColors=null;
  private static String MOL2IMG_SERVLETURL=null;

  /////////////////////////////////////////////////////////////////////////////
  public void doPost(HttpServletRequest request,HttpServletResponse response)
      throws IOException,ServletException
  {
    //SERVERPORT=request.getServerPort();
    SERVERNAME=request.getServerName();
    if (SERVERNAME.equals("localhost")) SERVERNAME=InetAddress.getLocalHost().getHostAddress();
    REMOTEHOST=request.getHeader("X-Forwarded-For"); // client (original)
    if (REMOTEHOST!=null)
    {
      String[] addrs=Pattern.compile(",").split(REMOTEHOST);
      if (addrs.length>0) REMOTEHOST=addrs[addrs.length-1];
    }
    else
    {
      REMOTEHOST=request.getRemoteAddr(); // client (may be proxy)
    }
    rb=ResourceBundle.getBundle("LocalStrings",request.getLocale());

    MultipartRequest mrequest=null;
    if (request.getMethod().equalsIgnoreCase("POST"))
    {
      try { mrequest=new MultipartRequest(request,UPLOADDIR,10*1024*1024,"ISO-8859-1",
                                    new DefaultFileRenamePolicy()); }
      catch (IOException lEx) {
        this.getServletContext().log("not a valid MultipartRequest",lEx); }
    }

    // main logic:
    ArrayList<String> cssincludes = new ArrayList<String>(Arrays.asList(CONTEXTPATH+"/css/biocomp.css"));
    ArrayList<String> jsincludes = new ArrayList<String>(Arrays.asList(CONTEXTPATH+"/js/biocomp.js",CONTEXTPATH+"/js/ddtip.js"));
    boolean ok=initialize(request,mrequest);
    if (mrequest!=null)	//method=POST, normal operation
    {
      if (!ok)
      {
        response.setContentType("text/html");
        out=response.getWriter();
        out.println(HtmUtils.HeaderHtm(APPNAME, jsincludes, cssincludes, JavaScript(), "", color1, request));
        out.println(HtmUtils.FooterHtm(errors,true));
        return;
      }
      else if (mrequest.getParameter("depict").equals("TRUE"))
      {
        response.setContentType("text/html");
        out=response.getWriter();
        out.println(HtmUtils.HeaderHtm(APPNAME, jsincludes, cssincludes, JavaScript(), "", color1, request));
        out.println(FormHtm(mrequest,response));
        Depict(mrequest,response);
        out.println(HtmUtils.OutputHtm(outputs));
        out.println(HtmUtils.FooterHtm(errors,true));
      }
    }
    else
    {
      if (request.getParameter("help")!=null)	// GET method, help=TRUE
      {
        response.setContentType("text/html");
        out=response.getWriter();
        out.println(HtmUtils.HeaderHtm(APPNAME, jsincludes, cssincludes, JavaScript(), "", color1, request));
        out.println(HelpHtm());
        out.println(HtmUtils.FooterHtm(errors,true));
      }
      else if (request.getParameter("test")!=null)	// GET method, test=TRUE
      {
        response.setContentType("text/plain");
        out=response.getWriter();
        HashMap<String,String> t = new HashMap<String,String>();
        t.put("JCHEM_LICENSE_OK",(LicenseManager.isLicensed(LicenseManager.JCHEM)?"True":"False"));
        t.put("JCHEM_MOLSEARCH_LICENSE_OK",(((new MolSearch()).isLicensed())?"True":"False"));
        out.print(HtmUtils.TestTxt(APPNAME,t));
      }
      else	// GET method, initial invocation of servlet w/ no params
      {
        response.setContentType("text/html");
        out=response.getWriter();
        out.println(HtmUtils.HeaderHtm(APPNAME, jsincludes, cssincludes, JavaScript(), "", color1, request));
        out.println(FormHtm(mrequest,response));
        out.println("<SCRIPT>go_init(window.document.mainform)</SCRIPT>");
        out.println(HtmUtils.FooterHtm(errors,true));
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  private boolean initialize(HttpServletRequest request,MultipartRequest mrequest)
      throws IOException,ServletException
  {
    SERVLETNAME=this.getServletName();

    MOL2IMG_SERVLETURL=(CONTEXTPATH+"/mol2img");

    outputs=new ArrayList<String>();
    errors=new ArrayList<String>();
    params=new HttpParams();
    mols=new ArrayList<Molecule>();
    sizes_h=new LinkedHashMap<String,Integer>();
    sizes_w=new LinkedHashMap<String,Integer>();

    String logo_htm="<TABLE CELLSPACING=5 CELLPADDING=5><TR><TD>";
    String imghtm=("<IMG BORDER=0 SRC=\""+CONTEXTPATH+"/images/biocomp_logo_only.gif\">");
    String tiphtm=(APPNAME+" web app from UNM Translational Informatics.");
    String href=("http://medicine.unm.edu/informatics/");
    logo_htm+=(HtmUtils.HtmTipper(imghtm,tiphtm,href,200,"white"));
    logo_htm+="</TD><TD>";
    imghtm=("<IMG BORDER=0 SRC=\""+CONTEXTPATH+"/images/chemaxon_powered_100px.png\">");
    tiphtm=("JChem from ChemAxon Ltd.");
    href=("http://www.chemaxon.com");
    logo_htm+=(HtmUtils.HtmTipper(imghtm,tiphtm,href,200,"white"));
    logo_htm+="</TD></TR></TABLE>";
    errors.add(logo_htm);

    //booleans:
    sizes_h.put("xs",96); sizes_w.put("xs",96);
    sizes_h.put("s",160); sizes_w.put("s",160);
    sizes_h.put("m",180); sizes_w.put("m",260);
    sizes_h.put("l",280); sizes_w.put("l",380);
    sizes_h.put("xl",480); sizes_w.put("xl",640);

    //Create webapp-specific log dir if necessary:
    File dout=new File(LOGDIR);
    if (!dout.exists())
    {
      boolean ok=dout.mkdir();
      System.err.println("LOGDIR creation "+(ok?"succeeded":"failed")+": "+LOGDIR);
      if (!ok)
      {
        errors.add("ERROR: could not create LOGDIR: "+LOGDIR);
        return false;
      }
    }

    String logpath=LOGDIR+"/"+SERVLETNAME+".log";
    logfile=new File(logpath);
    if (!logfile.exists())
    {
      try {
        logfile.createNewFile();
      }
      catch (IOException e)
      {
        errors.add("ERROR: Cannot create log file:"+e.getMessage());
        return false;
      }
      logfile.setWritable(true,true);
      PrintWriter out_log=new PrintWriter(logfile);
      out_log.println("date\tip\tN"); 
      out_log.flush();
      out_log.close();
    }
    if (!logfile.canWrite())
    {
      errors.add("ERROR: Log file not writable.");
      return false;
    }
    BufferedReader buff=new BufferedReader(new FileReader(logfile));
    if (buff==null)
    {
      errors.add("ERROR: Cannot open log file.");
      return false;
    }
    int n_lines=0;
    String line=null;
    String startdate=null;
    while ((line=buff.readLine())!=null)
    {
      ++n_lines;
      String[] fields=Pattern.compile("\\t").split(line);
      if (n_lines==2) startdate=fields[0];
    }
    buff.close(); //Else can result in error: "Too many open files"
    Calendar calendar=Calendar.getInstance();
    if (n_lines>2)
    {
      calendar.set(Integer.parseInt(startdate.substring(0,4)),
               Integer.parseInt(startdate.substring(4,6))-1,
               Integer.parseInt(startdate.substring(6,8)),
               Integer.parseInt(startdate.substring(8,10)),
               Integer.parseInt(startdate.substring(10,12)),0);

      DateFormat df=DateFormat.getDateInstance(DateFormat.FULL,Locale.US);
      errors.add("since "+df.format(calendar.getTime())+", times used: "+(n_lines-1));
    }

    calendar.setTime(new Date());
    DATESTR=String.format("%04d%02d%02d%02d%02d",
      calendar.get(Calendar.YEAR),
      calendar.get(Calendar.MONTH)+1,
      calendar.get(Calendar.DAY_OF_MONTH),
      calendar.get(Calendar.HOUR_OF_DAY),
      calendar.get(Calendar.MINUTE));

    LicenseManager.refresh();

    if (mrequest==null) return false;

    /// Stuff for a run:

    for (Enumeration e=mrequest.getParameterNames(); e.hasMoreElements(); )
    {
      String key=(String)e.nextElement();
      if (mrequest.getParameter(key)!=null)
        params.setVal(key,mrequest.getParameter(key));
    }

    if (params.isChecked("verbose"))
    {
      //errors.add("JChem version: "+chemaxon.jchem.version.VersionInfo.getVersion());
      errors.add("JChem version: "+com.chemaxon.version.VersionInfo.getVersion());
      errors.add("server: "+CONTEXT.getServerInfo()+" [API:"+CONTEXT.getMajorVersion()+"."+CONTEXT.getMinorVersion()+"]");
      errors.add("ServletName: "+this.getServletName());
    }

    String fname="infile";
    File ifile=mrequest.getFile(fname);
    String intxt=params.getVal("intxt").replaceFirst("[\\s]+$","");
    byte[] inbytes=new byte[1024];
    if (ifile!=null)
    {
      FileInputStream fis=new FileInputStream(ifile);
      int asize=inbytes.length;
      int size=0;
      int b;
      while ((b=fis.read())>=0)
      {
        if (size+1>asize)
        {
          asize*=2;
          byte[] tmp=new byte[asize];
          System.arraycopy(inbytes,0,tmp,0,size);
          inbytes=tmp;
        }
        inbytes[size]=(byte)b;
        ++size; 
      }
      byte[] tmp=new byte[size];
      System.arraycopy(inbytes,0,tmp,0,size);
      inbytes=tmp;
    }
    else if (intxt.length()>0)
    {
      if (params.getVal("molfmt").equals("cdx"))
      {
        inbytes=Base64Decoder.decodeToBytes(intxt);
      }
      else
      {
        inbytes=intxt.getBytes("utf-8");
      }
    }
    else
    {
      errors.add("No input data.");
      return false;
    }

    MolImporter molReader=null;
    if (params.getVal("molfmt").equals("automatic"))
    {
      String orig_fname=mrequest.getOriginalFileName(fname);
      String molfmt_auto=MFileFormatUtil.getMostLikelyMolFormat(orig_fname);
      if (orig_fname!=null && molfmt_auto!=null)
      {
        molReader=new MolImporter(new ByteArrayInputStream(inbytes),molfmt_auto);
      }
      else
      {
        molReader=new MolImporter(new ByteArrayInputStream(inbytes));
      }
    }
    else
    {
      String ifmt=params.getVal("molfmt");
      if (ifmt.startsWith("smiles") && params.isChecked("showprops"))
      {
        ifmt+=":";
        ArrayList<String> tags=SmiPropTags(inbytes);
        for (int i=0;i<tags.size();++i)
        {
          if (i>0) ifmt+=(",");
          ifmt+=("f"+tags.get(i));
        }
      }
      molReader=new MolImporter(new ByteArrayInputStream(inbytes),ifmt);
    }
    String fmt=molReader.getFormat();
    params.setVal("molfmt_auto",fmt);

    if (params.getVal("molfmt").equals("mrv") ||
             params.getVal("molfmt_auto").equals("mrv"))
    {
      // atomcolors idxs correspond to the mrvSetSeq values, assigned by MolAtom::setSetSeq().
      atomColors= new ArrayList<Color>();
      MolImporter molReader2=null;
      if (ifile!=null)
        molReader2 = new MolImporter(new FileInputStream(ifile));
      else
        molReader2=new MolImporter(new ByteArrayInputStream(inbytes),"mrv");
      MDocument mdoc=molReader2.nextDoc();
      molReader2.close();
      for (int i=0;i<10;++i)
      {
        Color c=mdoc.getAtomSetColor(i);
        if (c==null) break;
        atomColors.add(c);
        if (params.isChecked("verbose"))
          errors.add(String.format("atomcolor[%d]: %02X%02X%02X",i,c.getRed(),c.getGreen(),c.getBlue()));
      }
    }

    if (ifile!=null) ifile.delete();

    MFileFormat mffmt=MFileFormatUtil.getFormat(fmt);

    if (params.isChecked("file2txt"))
    {
      if (mffmt==MFileFormat.CDX) //binary
      {
        intxt=Base64Encoder.encode(inbytes);
        if (params.getVal("molfmt").equals("automatic"))
          params.setVal("molfmt","cdx");
      }
      else
      {
        intxt=new String(inbytes,"utf-8");
      }
      params.setVal("intxt",intxt);
    }

    Molecule m;
    int n2d=0;
    int n_failed=0;
    while (true)
    {
      try { m=molReader.read(); }
      catch (MolFormatException e)
      {
        errors.add("ERROR: MolImporter failed: "+e.getMessage());
        ++n_failed;
        continue;
      }
      if (m==null) break;
      if (m.getDim()==2) ++n2d;
      if (params.getVal("smarts").length()>0)
      {
        m.aromatize(MoleculeGraph.AROM_GENERAL); // aromatize so smarts work correctly.
        if (params.isChecked("align2smarts")&&(!params.isChecked("use2d")||m.getDim()<2))
          m.clean(2,null,null);
      }
      mols.add(m);
    }
    molReader.close();

    if (params.isChecked("parts2mols")) 
    {
      mols=Parts2mols(mols);
    }

    if (params.getVal("smarts").length()>0)
    {
      if (!(new MolSearch()).isLicensed())
      {
        errors.add("Warning: license error; smarts matching disabled.");
        molsearch=null;
      }
      else
      {
        MolHandler smartsReader=new MolHandler();
        molsearch=new MolSearch();
        smartsReader.setQueryMode(true);
        try {
          smartsReader.setMolecule(params.getVal("smarts"));
        }
        catch (MolFormatException e) {
          errors.add("ERROR: "+e.getMessage());
          errors.add("ERROR: smarts bad; matching disabled.");
          molsearch=null;
        }
        if (molsearch!=null)
        {
          try {
            molsearch.setQuery(smartsReader.getMolecule());
          }
          catch (LicenseException e) {
            errors.add("ERROR: "+e.getMessage());
            molsearch=null;
          }
        }
      }
    }
    if (params.isChecked("verbose"))
    {
      String desc=MFileFormatUtil.getFormat(molReader.getFormat()).getDescription();
      errors.add("input format:  "+molReader.getFormat()+" ("+desc+")");
      errors.add("mols read:  "+mols.size());
    }
    if (n_failed>0) errors.add("ERRORS (unable to read mol): "+n_failed);
    if (params.isChecked("use2d"))
      errors.add(String.format("%d / %d mols include 2D.",n2d,mols.size()));
    return true;
  }
  /////////////////////////////////////////////////////////////////////////////
  private static String FormHtm(MultipartRequest mrequest,HttpServletResponse response)
      throws IOException
  {
    String molfmt_menu="<SELECT NAME=\"molfmt\">\n";
    molfmt_menu+=("<OPTION VALUE=\"automatic\">automatic\n");
    for (String fmt: MFileFormatUtil.getMolfileFormats())
    {
      String desc=MFileFormatUtil.getFormat(fmt).getDescription();
      molfmt_menu+=("<OPTION VALUE=\""+fmt+"\">"+desc+"\n");
    }
    molfmt_menu+=("</SELECT>\n");
    molfmt_menu=molfmt_menu.replace(params.getVal("molfmt")+"\">",params.getVal("molfmt")+"\" SELECTED>\n");

    String size_menu="<SELECT NAME=\"size\">\n";
    for (String key:sizes_h.keySet())
    {
      size_menu+=("<OPTION VALUE=\""+key+"\">"+key+" - "+sizes_h.get(key)+"x"+sizes_w.get(key)+"\n");
    }
    size_menu+="</SELECT>\n";
    size_menu=size_menu.replace("\""+params.getVal("size")+"\">","\""+params.getVal("size")+"\" SELECTED>\n");

    String mode_menu="<SELECT NAME=\"mode\">\n";
    mode_menu+="<OPTION VALUE=\"bow\">BOW\n";
    mode_menu+="<OPTION VALUE=\"cob\">COB\n";
    mode_menu+="<OPTION VALUE=\"cow\">COW\n";
    mode_menu+="</SELECT>\n";
    mode_menu=mode_menu.replace("\""+params.getVal("mode")+"\">","\""+params.getVal("mode")+"\" SELECTED>\n");

    String ncols_menu="<SELECT NAME=\"ncols\">";
    ncols_menu+=("<OPTION VALUE=\"auto\">auto");
    for (int i=1;i<11;++i) ncols_menu+=("<OPTION VALUE=\""+i+"\">"+i);
    ncols_menu+="</SELECT>\n";
    ncols_menu=ncols_menu.replace("\""+params.getVal("ncols")+"\">","\""+params.getVal("ncols")+"\" SELECTED>\n");

    String imgfmt_png="";
    String imgfmt_jpeg="";
    if (params.getVal("imgfmt").equals("jpeg")) imgfmt_jpeg="CHECKED";
    else imgfmt_png="CHECKED";

    String htm=
    ("<FORM NAME=\"mainform\" METHOD=POST")
    +(" ACTION=\""+response.encodeURL(SERVLETNAME)+"\"")
    +(" ENCTYPE=\"multipart/form-data\">\n")
    +("<INPUT TYPE=HIDDEN NAME=\"depict\">\n")
    +("<TABLE WIDTH=\"100%\"><TR><TD><H1>"+APPNAME+"</H1></TD>\n")
    +("<TD ALIGN=RIGHT>\n")
    +("<BUTTON TYPE=BUTTON onClick=\"void window.open('"+response.encodeURL(SERVLETNAME)+"?help=TRUE','helpwin','width=600,height=400,scrollbars=1,resizable=1')\"><B>Help</B></BUTTON>\n")
    +("<BUTTON TYPE=BUTTON onClick=\"go_demo(this.form)\"><B>Demo</B></BUTTON>\n")
    +("<BUTTON TYPE=BUTTON onClick=\"window.location.replace('"+response.encodeURL(SERVLETNAME)+"')\"><B>Reset</B></BUTTON>\n")
    +("</TD></TR></TABLE>\n")
    +("<HR>\n")
    +("<TABLE WIDTH=\"100%\" CELLPADDING=5 CELLSPACING=5>\n")
    +("<TR BGCOLOR=\"#CCCCCC\"><TD VALIGN=TOP>\n")
    +("format:"+molfmt_menu)
    +("<INPUT TYPE=CHECKBOX NAME=\"file2txt\" VALUE=\"CHECKED\"")
    +(" "+params.getVal("file2txt")+">file2txt<BR>\n")
    +("upload: <INPUT TYPE=\"FILE\" NAME=\"infile\"> ...or paste:")
    +("<BR><TEXTAREA NAME=\"intxt\" WRAP=OFF ROWS=12 COLS=60>")
    +(params.getVal("intxt"))
    +("</TEXTAREA>\n")
    +("</TD>\n")
    +("<TD VALIGN=TOP>\n")
    +("<B>options:</B><BR>\n")
    +("<TABLE WIDTH=100%><TR><TD VALIGN=TOP>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"showh\" VALUE=\"CHECKED\"")
    +(params.getVal("showh")+">show Hs<BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"lonepairs\" VALUE=\"CHECKED\"")
    +(params.getVal("lonepairs")+">lonepairs<BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"use2d\" VALUE=\"CHECKED\"")
    +(params.getVal("use2d")+">use2d<BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"transparent\" VALUE=\"CHECKED\"")
    +(" "+params.getVal("transparent")+">transparent<BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"zoomable\" VALUE=\"CHECKED\"")
    +(" "+params.getVal("zoomable")+">zoomable<BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"verbose\" VALUE=\"CHECKED\"")
    +(" "+params.getVal("verbose")+">verbose<BR>\n")
    +("</TD><TD VALIGN=TOP>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"showprops\" VALUE=\"CHECKED\"")
    +(" "+params.getVal("showprops")+">show properties<BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"showmaps\" VALUE=\"CHECKED\"")
    +(" "+params.getVal("showmaps")+">show atom maps<BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"smilesmatch\" VALUE=\"CHECKED\"")
    +(" "+params.getVal("smilesmatch")+">smilesmatch<BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"parts2mols\" VALUE=\"CHECKED\"")
    +(" "+params.getVal("parts2mols")+">parts2mols<BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"showarom\" VALUE=\"CHECKED\"")
    +(" "+params.getVal("showarom")+">showarom<BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"align2smarts\" VALUE=\"CHECKED\"")
    +(" "+params.getVal("align2smarts")+">align2smarts<BR>\n")
    +("</TD></TR>\n")
    +("</TABLE>\n")
    +("smarts: <INPUT TYPE=\"text\" NAME=\"smarts\" SIZE=\"48\"")
    +(" VALUE=\""+params.getVal("smarts")+"\">\n")
    +("<BR>\n")
    +("<B>output:</B><BR>\n")
    +("size:"+size_menu)
    +("ncols:"+ncols_menu)
    +("mode:"+mode_menu+"<BR>\n")
    +("format: <INPUT TYPE=RADIO NAME=\"imgfmt\" VALUE=\"png\" "+imgfmt_png+">png")
    +("&nbsp;<INPUT TYPE=RADIO NAME=\"imgfmt\" VALUE=\"jpeg\" "+imgfmt_jpeg+">jpeg<BR>\n")
    +("</TD></TR></TABLE>\n")
    +("<P>\n")
    +("<CENTER>\n")
    +("<BUTTON TYPE=BUTTON onClick=\"go_depict(this.form)\"><B>Go "+APPNAME+"</B></BUTTON>\n")
    +("</CENTER>\n")
    +("</FORM>\n");
    return htm;
  }
  /////////////////////////////////////////////////////////////////////////////
  private static ArrayList<Molecule> Parts2mols(ArrayList<Molecule> mols)
  {
    ArrayList<Molecule> partmols=new ArrayList<Molecule>();
    for (Molecule mol:mols)
    {
      Molecule[] partmols_this=mol.convertToFrags();
      int i_part=0;
      for (Molecule partmol:partmols_this)
      {
        ++i_part;
        partmol.setName(mol.getName()+" ["+i_part+"]");
        partmols.add(partmol);
      } 
    }
    return partmols;
  }
  /////////////////////////////////////////////////////////////////////////////
  private static void Depict(MultipartRequest mrequest,HttpServletResponse response)
      throws IOException,ServletException
  {
    int n_mols=0;
    int w=sizes_w.get(params.getVal("size"));
    int h=sizes_h.get(params.getVal("size"));
    int n_cols=1;
    if (params.getVal("ncols").equals("auto"))
      n_cols=(int)(900.0/(float)w);
    else
      n_cols=Integer.parseInt(params.getVal("ncols"));
    String depictopts=("mode="+params.getVal("mode"));
    if (params.getVal("imgfmt").equals("jpeg")) depictopts+=("&imgfmt=jpeg");
    else depictopts+=("&imgfmt=png");
    if (params.isChecked("showarom")) depictopts+=("&arom_gen=true");
    else depictopts+=("&kekule=true");
    if (params.isChecked("showh")) depictopts+=("&showh=true");
    if (params.isChecked("showmaps")) depictopts+=("&showmaps=true");
    if (params.isChecked("lonepairs")) depictopts+=("&lonepairs=true");
    if (params.isChecked("use2d")) depictopts+=("&use2d=true");
    if (params.isChecked("transparent")) depictopts+=("&transparent=true");
    if (params.getVal("smarts").length()>0 && molsearch!=null)
    {
      depictopts+=("&smartscode="+URLEncoder.encode(params.getVal("smarts"),"UTF-8"));
      if (params.isChecked("align2smarts"))
      {
        int n_aligned=molalign_utils.AlignToSmarts(mols,params.getVal("smarts"));
        if (params.isChecked("verbose"))
          errors.add("aligned: "+n_aligned);
      }
      if (params.isChecked("smilesmatch")) depictopts+=("&smilesmatch=true");
    }
    if (!params.isChecked("showprops"))  // properties means both SD data and query properties.
      depictopts+=("&clearqprops=TRUE");

    depictopts+=("&maxscale=0"); //needed for JChem

    String smiles=null;
    String mdlcode=null;
    String mrvcode=null;
    String molname=null;
    String prophtm=("<OL>\n");
    String tablehtm=("<CENTER><TABLE BORDER>\n");
    for (Molecule mol:mols)
    {
      if (n_mols%n_cols==0) tablehtm+="<TR>\n";
      if (params.getVal("molfmt").equals("mrv") ||
               params.getVal("molfmt_auto").equals("mrv"))
      {
        try {
          if (params.isChecked("showh")) 
            mrvcode=mol.exportToFormat("base64:gzip:mrv:+H");
          else
            mrvcode=mol.exportToFormat("base64:gzip:mrv");
        }
        catch (MolExportException e) {
          errors.add(e.getMessage());
          return;
        }
      }
      else if (params.isChecked("use2d") && mol.getDim()==2)
      {
        try {
          if (params.isChecked("showh")) 
            mdlcode=mol.exportToFormat("base64:gzip:mol:+H");
          else
            mdlcode=mol.exportToFormat("base64:gzip:mol:-H");
        }
        catch (MolExportException e) {
          errors.add(e.getMessage());
          return;
        }
      }
      else if (params.getVal("molfmt").equals("smarts") ||
               params.getVal("molfmt_auto").equals("smarts"))
      {
        // Note: Cannot export smarts query-mol as smiles;
        // however, can depict a query-mol.
        smiles=mol.exportToFormat("smarts");
      }
      else
      {
        try { smiles=mol.exportToFormat("smiles:u"); }
        catch (MolExportException e) {
          if (params.isChecked("verbose")) errors.add(e.getMessage());
          try { smiles=mol.exportToFormat("cxsmiles:u"); }
          catch (MolExportException e2) {
            if (params.isChecked("verbose")) errors.add(e2.getMessage());
            try { smiles=mol.exportToFormat("smarts:u"); }
            catch (MolExportException e3) {
              if (params.isChecked("verbose")) errors.add(e3.getMessage());
              try { smiles=mol.exportToFormat("cxsmarts:u"); }
              catch (MolExportException e4) {
                errors.add(e4.getMessage());
                smiles="";
              }
            }
          }
        }
      }

      String prophtm_this="";
      if (params.isChecked("showprops"))
      {
        prophtm_this+=("<UL>\n");
        for  (String key:mol.properties().getKeys())
        {
          prophtm_this+=("<LI>"+key+":"+mol.getProperty(key)+"\n");
        }
        prophtm_this+=("</UL>\n");
        prophtm+=("<LI>"+mol.getName()+"\n"+prophtm_this);
      }

      String imhtm="";
      if (mdlcode!=null)
      {
        imhtm=HtmUtils.Mdlcode2ImgHtm(mdlcode,depictopts,h,w,
                          MOL2IMG_SERVLETURL,
                          params.getVal("zoomable").equals("CHECKED"),4,
                          "go_zoom_mdl2img");
      }
      else if (mrvcode!=null)
      {
        imhtm=HtmUtils.Mrvcode2ImgHtm(mrvcode,atomColors,depictopts,h,w,
                          MOL2IMG_SERVLETURL,
                          params.getVal("zoomable").equals("CHECKED"),4,
                          "go_zoom_mrv2img");
      }
      else
      {
        imhtm=HtmUtils.Smi2ImgHtm(smiles,depictopts,h,w,
                          MOL2IMG_SERVLETURL,
                          params.getVal("zoomable").equals("CHECKED"),4,
                          "go_zoom_smi2img");
      }
      tablehtm+=("<TD BGCOLOR=\"white\" ALIGN=\"CENTER\">"+imhtm+"<BR>");
      tablehtm+=("<TT>");
      tablehtm+=(mol.getName());
      if (params.isChecked("showprops"))
        tablehtm+=("<BR>\n"+prophtm_this);
      tablehtm+=("</TT></TD>\n");
      ++n_mols;
      if (params.isChecked("verbose"))
      {
        if (!mol.getName().isEmpty()) errors.add(""+n_mols+". "+mol.getName());
      }

      if (params.getVal("smarts").length()>0 && params.isChecked("verbose") && molsearch!=null)
      {
        molsearch.setTarget(mol);
        int[][] matchs=null;
        try {
          matchs=molsearch.findAll();
        }
        catch (SearchException e) {
          errors.add(e.getMessage());
        }
        int n_matchs=((matchs!=null)?matchs.length:0);
        errors.add("smarts matches: "+n_matchs);
      }

      if (n_mols%n_cols==0) tablehtm+="</TR>\n";
      if (n_mols==N_MAX) break;
    }
    if (n_mols%n_cols>0)
    {
      if (n_mols>n_cols)
      {
        for (int i=n_mols%n_cols;i<n_cols;++i)
        { tablehtm+=("<TD ALIGN=CENTER>~</TD>\n"); }
      }
      tablehtm+=("</TR>");
    }
    tablehtm+=("</TABLE></CENTER>");
    outputs.add(tablehtm);
    outputs.add(prophtm+"</OL>\n");
    errors.add("n_mols: "+n_mols);

    PrintWriter out_log=new PrintWriter(
      new BufferedWriter(new FileWriter(logfile,true)));
    out_log.printf("%s\t%s\t%d\n",DATESTR,REMOTEHOST,n_mols); 
    out_log.close();
  }
  /////////////////////////////////////////////////////////////////////////////
  private static ArrayList<String> SmiPropTags(byte[] inbytes)
    throws IOException,UnsupportedEncodingException
  {
    ArrayList<String> tags=new ArrayList<String>();
    int nfields=0;
    String intxt=new String(inbytes,"utf-8");
    BufferedReader buff=new BufferedReader(new StringReader(intxt));
    String line=buff.readLine();
    if (line!=null)
    {
      String[] fields=Pattern.compile("\\t").split(line);
      nfields=fields.length-1;
      if (fields[0].startsWith("#"))
      {
        fields[0]=fields[0].replace("^#+","");
        for (int i=1;i<=nfields;++i)
          tags.add(fields[i]);
      }
      else
      {
        for (int i=1;i<=nfields;++i)
          tags.add(String.format("field%02d",i));
      }
    }
    buff.close();
    return tags;
  }
  /////////////////////////////////////////////////////////////////////////////
  private static String JavaScript()
  {
    return(
"function go_init(form)"+
"{\n"+
"  form.file2txt.checked=true;\n"+
"  form.smarts.value='';\n"+
"  form.intxt.value='';\n"+
"  var i;\n"+
"  for (i=0;i<form.molfmt.length;++i)\n"+
"    if (form.molfmt.options[i].value=='automatic')\n"+
"      form.molfmt.options[i].selected=true;\n"+
"  for (i=0;i<form.size.length;++i)\n"+
"    if (form.size.options[i].value=='m')\n"+
"      form.size.options[i].selected=true;\n"+
"  for (i=0;i<form.imgfmt.length;++i)\n"+ //radio
"    if (form.imgfmt[i].value=='png')\n"+
"      form.imgfmt[i].checked=true;\n"+
"  for (i=0;i<form.mode.length;++i)\n"+
"    if (form.mode.options[i].value=='cow')\n"+
"      form.mode.options[i].selected=true;\n"+
"  form.showh.checked=false;\n"+
"  form.showarom.checked=false;\n"+
"  form.lonepairs.checked=false;\n"+
"  form.use2d.checked=true;\n"+
"  form.zoomable.checked=true;\n"+
"  form.showprops.checked=false;\n"+
"  form.transparent.checked=false;\n"+
"  form.verbose.checked=false;\n"+
"}\n"+
"function checkform(form)\n"+
"{\n"+
"  if (!form.intxt.value && !form.infile.value) {\n"+
"    alert('ERROR: No input specified');\n"+
"    return 0;\n"+
"  }\n"+
"  return 1;\n"+
"}\n"+
"function go_demo(form)\n"+
"{\n"+
"  go_init(form);"+
"  form.intxt.value='CN1C(=O)N(C)C(=O)C(N(C)C=N2)=C12 caffeine\\n';\n"+
"  form.intxt.value+='COc1cc2c(ccnc2cc1)C(O)C4CC(CC3)C(C=C)CN34 quinine\\n';\n"+
"  form.intxt.value+='CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O benzylpenicillin\\n';\n"+
"  form.intxt.value+='CCC(=C(c1ccc(OCCN(C)C)cc1)c1ccccc1)c1ccccc1 Tamoxifen\\n';\n"+
"  form.intxt.value+='CNCCC(c1ccccc1)Oc2ccc(cc2)C(F)(F)F.Cl Prozac\\n';\n"+
"  form.intxt.value+='NC(C)Cc1ccccc1 adderall\\n';\n"+
"  form.intxt.value+='CNC(=C[N+](=O)[O-])NCCSCC1=CC=C(O1)CN(C)C.Cl Zantac\\n';\n"+
"  form.intxt.value+='Oc2cc(cc1OC(C3CCC(=CC3c12)C)(C)C)CCCCC THC\\n';\n"+
"  form.intxt.value+='CN1c2ccc(cc2C(=NCC1=O)c3ccccc3)Cl Valium\\n';\n"+
"  form.depict.value='TRUE'\n"+
"  form.submit()\n"+
"}\n"+
"function go_depict(form)\n"+
"{\n"+
"  if (!checkform(form)) return;\n"+
"  form.depict.value='TRUE'\n"+
"  form.submit()\n"+
"}\n");
  }
  /////////////////////////////////////////////////////////////////////////////
  private static String HelpHtm()
  {
    return (
    "<B>"+APPNAME+" help</B><P>\n"+
    "This web app consists of (1) a Java servlet using JChem\n"+
    "for the user interface, and (2) a Java servlet using JChem which\n"+
    "generates inline graphics, PNG or JPEG.\n"+
    "Molecules with associated 2D coordinates should be\n"+
    "depicted using those coordinates.\n"+
    "The list of molecule formats is automatically generated by JChem.\n"+
    "<P>\n"+
    "Note that show properties pertains to \n"+
    "(1) smiles files with multicolumn data, with or without an\n"+
    "initial #-prefixed header/comment line, and (2) SDF, RDF,\n"+
    "PDB or other files which contain associated data.\n"+
    "<P>\n"+
    "The smilesmatch option causes the input smarts to be parsed\n"+
    "as a smiles.  One result is that explicitly denoted hydrogens (e.g. \"[nH]\") are\n"+
    "not interpreted as query specifications.\n"+
    "<P>\n"+
    "configured with molecule limit N_MAX = "+N_MAX+"\n"+
    "<P>\n"+
    "Thanks to ChemAxon for the use of JChem in this application.\n"+
    "<P>\n"+
    "author/support: Jeremy Yang\n"
    );
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Read servlet parameters (from web.xml).
  */
  public void init(ServletConfig conf) throws ServletException
  {
    super.init(conf);
    CONTEXT=getServletContext();	// inherited method
    CONTEXTPATH=CONTEXT.getContextPath();
    try { APPNAME=conf.getInitParameter("APPNAME"); }
    catch (Exception e) { APPNAME=this.getServletName(); }
    UPLOADDIR=conf.getInitParameter("UPLOADDIR");
    if (UPLOADDIR==null)
      throw new ServletException("Please supply UPLOADDIR parameter");
    LOGDIR=conf.getInitParameter("LOGDIR")+CONTEXTPATH;
    if (LOGDIR==null) LOGDIR="/tmp"+CONTEXTPATH+"_logs";
    try { N_MAX=Integer.parseInt(conf.getInitParameter("N_MAX")); }
    catch (Exception e) { N_MAX=100; }
  }
  /////////////////////////////////////////////////////////////////////////////
  public void doGet(HttpServletRequest request,HttpServletResponse response)
      throws IOException, ServletException
  {
    doPost(request,response);
  }
}
