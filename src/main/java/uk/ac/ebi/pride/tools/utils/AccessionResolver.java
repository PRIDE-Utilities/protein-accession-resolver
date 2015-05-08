package uk.ac.ebi.pride.tools.utils;

import org.apache.log4j.Logger;

import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * User: rcote
 * Date: 08-Nov-2007
 * Time: 12:42:12
 */
public class AccessionResolver {

    private static Logger logger = Logger.getLogger(AccessionResolver.class);

    //GI patterm
    private static final Pattern GI_PATTERN = Pattern.compile("^(GI)? ?\\|? ?([0-9]+)(\\|)?.*");
    private static final Pattern SP_PATTERN = Pattern.compile("^(SP|TR|TRM) ?\\|? ?([A-Z][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9])(-[0-9]+)(\\|)?.*");
    private static final Pattern SIMPLE_SP = Pattern.compile("([A-Z][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9])");
    private static final Pattern MIDDLE_PIPE = Pattern.compile(".*\\|(.*)\\|.*");
    private static final Pattern BAD_IPI_PATTERN = Pattern.compile("IPI[0-9]*([OPQ][0-9A-Z]*-?[0-9A-Z]*)");
    private static final Pattern BAD_PATTERN_1 = Pattern.compile("[0-9]+/[0-9]+");
    private static final Pattern UNIREF_PATTERN = Pattern.compile("UNIREF[0-9]*_?([A-Z][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9])");
    private static final long MIN_VALID_GI = 1000;

    //minimum length of accession to send to PICR
    private static final int MIN_VALID_AC_LENGTH = 5;

    //accession and version that are parsed out
    private String accession;
    private String version;
    //sanity check
    boolean isValidAccession = true;

    private String accessionToMap;
    private String versionToMap;
    private String database;
    private boolean isHybridDatabase;

    /**
     * the constructor without database name
     */
    public AccessionResolver(String accessionToMap, String versionToMap) {
        this(accessionToMap, versionToMap, "");
    }

    /**
     * the constructor of uk.ac.ebi.pride.tools.utils.AccessionResolver will try and extract a
     * valid accession and version from whatever was submitted to PRIDE.
     */
    public AccessionResolver(String accessionToMap, String versionToMap, String database) {
        this(accessionToMap, versionToMap, database, false);
    }

    /**
     * the constructor of uk.ac.ebi.pride.tools.utils.AccessionResolver will try and extract a
     * valid accession and version from whatever was submitted to PRIDE.
     */
    public AccessionResolver(String accessionToMap, String versionToMap, String database, boolean isHybridDatabase) {
        this.accessionToMap = accessionToMap;
        this.versionToMap = versionToMap;
        this.database = database;
        this.isHybridDatabase = isHybridDatabase;

        if (this.accessionToMap == null || this.database == null) {
            throw new IllegalArgumentException("Please provide a valid, non-null, accession and database name to the accession resolver");
        }

        fixMalformedAccessionAndVersion();

        if (!isValidAccession) {
            StringBuilder message = new StringBuilder();
            message.append("INVALID").append('\t')
                    .append(accessionToMap).append('\t')
                    .append(versionToMap).append('\t')
                    .append(database).append('\t')
                    .append("PARSED ->").append(accession).append("<-");
            logger.debug(message.toString());
        }
    }

    public String getAccession() {
        return accession;
    }

    public String getVersion() {
        return version;
    }

    public boolean isValidAccession() {
        return isValidAccession;
    }

    private void fixMalformedAccessionAndVersion() {

//submitted accessions in pride can be generally f'ed up. clean up before mapping
//        JQ1489 thymosin beta-4 - African clawed frog
//        JW0049 protein-tyrosine-phosphatase (EC 3.1.3.48)
//        O00264, NP_006658
//        P08670, NP_003371
//        P12105_1 [Segment 1 of 3] Collagen alpha 1(III) ch
//        O00154|BACH_HUMAN -- bad!
//        gi|1234234|BLA    -- good!
//        IPIUPI | IPIIPI | IPI[NXPOQ]

//reversed sequences should be ignored!
//      ENSDARP00000011638:reversed
//      gi|158275145|:reversed

//randomized sequences should be ignored!
//      ###RND###Q09472|EP300_HUMAN.

//two conflicting usecases for versioning!
//        IPI0001234..5     -- need to deal with version that is .5
//        Tb10.6k15.3810    -- need to recognize that there is no version

//for TAIR, AT3G17770.1 is the accession, there is no version


        //in the case of hybrid databases, you can have decoy and valid sequences mixed together.
        //if this is the case, skip this section as you can still map valid accessions from this database
        if (!isHybridDatabase) {

            //first sanity check on reversed databases
            if (database.toLowerCase().contains("reverse") && !database.toLowerCase().contains("_concat")) {
                isValidAccession = false;
                return;
            }
            if (database.toLowerCase().contains("_rev") && !database.toLowerCase().contains("_concat")) {
                isValidAccession = false;
                return;
            }
            if (database.toLowerCase().contains("decoy")) {
                isValidAccession = false;
                return;
            }

            if (database.toLowerCase().contains("dec_")) {
                isValidAccession = false;
                return;
            }

            if (database.toLowerCase().contains("_dec")) {
                isValidAccession = false;
                return;
            }

        }

        //start accession checks
        //if the accession contains "reverse" but not "reverse sense"
        if (accessionToMap.toLowerCase().contains("reverse") && !accessionToMap.toLowerCase().contains("reverse sense")) {
            isValidAccession = false;
            return;
        }

        if (accessionToMap.toLowerCase().contains("_rev")) {
            isValidAccession = false;
            return;
        }

        if (accessionToMap.toLowerCase().contains("rev_")) {
            isValidAccession = false;
            return;
        }

        if (accessionToMap.toLowerCase().contains("###rnd###")) {
            isValidAccession = false;
            return;
        }

        if (accessionToMap.toLowerCase().contains("###rev###")) {
            isValidAccession = false;
            return;
        }

        if (accessionToMap.toLowerCase().contains(".fasta")) {
            isValidAccession = false;
            return;
        }

        if (accessionToMap.toLowerCase().contains("jgi|aspni1")) {
            isValidAccession = false;
            return;
        }

        Matcher match = BAD_PATTERN_1.matcher(accessionToMap.toLowerCase());
        if (match.matches()) {
            isValidAccession = false;
            return;
        }

        //Q8C7X2-00-00-00 - why the hell not?!
        //Q9CXW4-00-00-00
        if (accessionToMap.endsWith("-00-00-00")) {
            accessionToMap = accessionToMap.substring(0, accessionToMap.indexOf("-00-00-00"));
        }

        //remove anything after the first space - i.e. use the first word as the ac
        int ndx = accessionToMap.indexOf(' ');
        if (ndx > 0) {

            boolean matchFound = false;

            match = SP_PATTERN.matcher(accessionToMap.toUpperCase());
            if (match.matches()) {
                accessionToMap = match.group(2);
                matchFound = true;
            }

            //check to make sure that it's a gi
            match = GI_PATTERN.matcher(accessionToMap.toUpperCase());
            if (match.matches()) {
                if (testGI(match.group(2))) {
                    accessionToMap = match.group(2);
                    matchFound = true;
                } else {
                    isValidAccession = false;
                    return;
                }
            }

            if (!matchFound) {
                //if it's not a GI/SP, trim the spaces and use first word
                accessionToMap = accessionToMap.substring(0, ndx).trim();
            }

        }
        //remove any possible commas
        ndx = accessionToMap.indexOf(',');
        if (ndx > 0) {
            accessionToMap = accessionToMap.substring(0, ndx).trim();
        }

        //check GIs or sp|P12345|BLA_HUMAN

        //remove any pipes but be careful of GIs!
        ndx = accessionToMap.indexOf('|');
        if (ndx > 0) {
            boolean matchFound = false;

            if (!matchFound) {
                //check to see if it's a SP
                match = SP_PATTERN.matcher(accessionToMap.toUpperCase());
                if (match.matches()) {
                    accessionToMap = match.group(2);
                    matchFound = true;
                }
            }

            if (!matchFound) {
                //check to see if it's a GI
                match = GI_PATTERN.matcher(accessionToMap.toUpperCase());
                if (match.matches()) {
                    if (testGI(match.group(2))) {
                        accessionToMap = match.group(2);
                        matchFound = true;
                    }
                }
            }

            if (!matchFound) {
                //try one last chance at matching patterns of the form
                //xxx|bla|yyy - generally we can find something useful for bla
                match = MIDDLE_PIPE.matcher(accessionToMap.toUpperCase());
                if (match.matches()) {

                    String tmpAc = match.group(1);
                    //is it a SP/TR?
                    match = SIMPLE_SP.matcher(tmpAc);
                    if (match.matches()) {
                        accessionToMap = tmpAc;
                        matchFound = true;
                    } else {
                        //is it a GI
                        if (testGI(tmpAc)) {
                            accessionToMap = tmpAc;
                            matchFound = true;
                        }
                    }
                }
            }

            if (!matchFound) {
                //O34528|yrvN
                //check to see if the first part of the pipe is a simple swissprot
                String tmpAc = accessionToMap.substring(0, ndx);
                match = SIMPLE_SP.matcher(tmpAc);
                if (match.matches()) {
                    accessionToMap = tmpAc;
                    matchFound = true;
                }
            }

            if (!matchFound) {
                //it's not a GI, it's not a SP
                //so trim the pipe and use whatever's left
                accessionToMap = accessionToMap.substring(ndx + 1);
            }

        }

        //check uniref
        match = UNIREF_PATTERN.matcher((accessionToMap.toUpperCase()));
        if (match.matches()) {
            accessionToMap = match.group(1);
        }

        //fix bad IPIIPI/IPIUPI/etc
        //        IPIUPI | IPIIPI | IPI[NXPOQ]
        HashSet<String> badIPIs = new HashSet<String>();
        badIPIs.add("IPII");
        badIPIs.add("IPIU");
        badIPIs.add("IPIN");
        badIPIs.add("IPIX");
        badIPIs.add("IPIP");
        badIPIs.add("IPIO");
        badIPIs.add("IPIQ");
        if (accessionToMap.startsWith("IPI")) {

            try {
                if (badIPIs.contains(accessionToMap.substring(0, 4))) {
                    accessionToMap = accessionToMap.substring(3);
                }
                match = BAD_IPI_PATTERN.matcher(accessionToMap);
                if (match.matches()) {
                    accessionToMap = match.group(1);
                }

            } catch (java.lang.StringIndexOutOfBoundsException saoob) {
                logger.error("Invalid IPI accession: " + accessionToMap);
                isValidAccession = false;
                return;
            }

        }

        //remove any quotes
        ndx = accessionToMap.indexOf("\"");
        if (ndx >= 0) {
            accessionToMap = accessionToMap.substring(ndx + 1);
        }
        ndx = accessionToMap.lastIndexOf("\"");
        if (ndx > 0) {
            accessionToMap = accessionToMap.substring(0, ndx);
        }

        //check for versions
        ndx = accessionToMap.indexOf('.');
        if (ndx > 0) {
            versionToMap = accessionToMap.substring(ndx + 1).trim();
            accessionToMap = accessionToMap.substring(0, ndx).trim();
        }

        //also fix weirdo versions - remove any dots and such
        //.1
        if (versionToMap != null && versionToMap.startsWith(".")) {
            ndx = versionToMap.lastIndexOf('.');
            versionToMap = versionToMap.substring(ndx + 1);
        }

        //try and check to see if we have valid versions
        //according to the above ckecks:
        //Tb10.6k15.3810 -> AC = Tb10
        //               -> Version = 6k15.3810
        try {

            //this will barf if the accession is bad.
            if (versionToMap != null) {
                new Integer(versionToMap);
            }

            //if we're dealing with a TAIR accession, this will restore
            //the incorrectly split accession
            if (database != null && database.toUpperCase().startsWith("TAIR") && versionToMap != null) {
                accessionToMap = accessionToMap + "." + versionToMap;
                versionToMap = null;
            }

            //if all is well, set proper versions
            accession = accessionToMap;
            version = versionToMap;

        } catch (NumberFormatException e) {

            accessionToMap = accessionToMap + "." + versionToMap;
            accession = accessionToMap;
            version = null;
            logger.info("Improper accession version detected, setting it to null. Accession to map is now: " + accessionToMap);

        }

        //final sanity check
        if (accessionToMap.length() < MIN_VALID_AC_LENGTH) {
            isValidAccession = false;
        }

    }

    private boolean testGI(String giString) {

        boolean retval = true;

        //if it is a GI, check to see that it's not less than 1000
        //this is an arbitrary value, but will catch most erroneous GI
        //patterns
        try {
            Long testGI = new Long(giString);
            if (testGI < MIN_VALID_GI) {
                retval = false;
            }
        } catch (NumberFormatException nfe) {
            //not a valid numerical GI
            retval = false;
        }

        return retval;

    }

    public static void main(String[] args) {

        AccessionResolver res = new AccessionResolver("npol_c_46_98 Neisseria polysaccharea ATCC 43768 [1688 - 654] (REVERSE SENSE) ORF", null, "human-crap-homd_TD.fasta.pro (no description)", true);
//        AccessionResolver res = new AccessionResolver("gi|108562577|ref|YP_626893.1|", null, "microhumancrapdecoycombo.fasta.pro", true);
        if (res.isValidAccession()) {
            System.out.println("res.getAccession() = " + res.getAccession());
            System.out.println("res.getVersion() = " + res.getVersion());
        }

    }

}
