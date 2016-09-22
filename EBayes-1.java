/**
 * EBayes.java
 * @author Fabio G. Cozman
 * Copyright 1998, 1999 Fabio Cozman, Universidade de Sao Paulo
 * fgcozman@usp.br, http://www.cs.cmu.edu/~fgcozman/home.html
 *
 * The EBayes distribution is free software; you can
 * redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, provided
 * that this notice and the name of the author appear in all copies.
 * If you're using the software, please notify fgcozman@usp.br so
 * that you can receive updates and patches. EBayes is distributed
 * "as is", in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with the EBayes distribution. If not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

import java.io.*;
import java.util.StringTokenizer;
import BayesianNetworks.*;
import BayesianInferences.*;

public class EBayes {
  BayesNet network;

  /**
   * Main method for EBayes.
   */
  public static void main(String args[]) {
    // Define network and input stream.
    EBayes ejb = new EBayes();
    DataInputStream dis = null;
    try {
        if (args.length > 0)
            dis = new DataInputStream(new FileInputStream(args[0]));
        else
            dis = new DataInputStream(System.in);
    } catch (IOException e) {
        System.exit(0);
    }

    // Announce with help message.
    help_message();

    // Read input and process queries.
    String line, command, name = null, value = null;
    StringTokenizer st;
    while(true) {
        try {
            System.out.print("Insert command character (l|o|t|u|i|m|e|x|q):\n>>");
            line = dis.readLine();
            st = new StringTokenizer(line);
            if (st.hasMoreElements()) {
                command = st.nextToken();
                if (st.hasMoreTokens()) {
                    name = st.nextToken();
                    if (st.hasMoreTokens())
                        value = st.nextToken();
                    else
                        value = null;
                } else
                    name = null;
                if (!ejb.process_command(command, name, value))
                    break;
            }
        } catch(IOException e) { }
    }
  }

  /**
   * Process the command inserted by the user.
   */
  boolean process_command(String command, String variable, String value) {
    int i;
    DiscreteVariable pv;
    Inference inf;
    Expectation exp;
    Explanation expl;
    DiscreteFunction result;
    double expectation;
    String[][] explanation;
    Class theNetworkClass;

    // Announce command.
    System.out.print("\tParsed: " + command + " ");
    if (variable != null) System.out.print(variable + " ");
    if (value != null) System.out.println(value);
    else               System.out.println();
    // Check that network is not null.
    if ((command.charAt(0) != 'l') && (command.charAt(0) != 'q') &&
        (command.charAt(0) != 'h') && (network == null)) {
        System.out.println("\tNetwork is not defined!");
        return(true);
    }
    // Perform operations.
    switch (command.charAt(0)) {
        case 'q': return(false);
        case 'h': help_message(); break;
        case 'l':
            try {
                theNetworkClass = Class.forName(variable);
                System.out.println("\tLoaded " + theNetworkClass);
                network = (BayesNet)(theNetworkClass.newInstance());
            } catch (ClassNotFoundException e) {
                System.out.println("\tDefining class was not found!");
            } catch (IllegalAccessException e) {
                System.out.println("\tNot possible to access defining class!");
            } catch (InstantiationException e) {
                System.out.println("\tDefining class cannot be instantiated!");
            }
            break;
        case 'o':
            if (variable == null) {
                System.out.println("\tNull variable name!");
                return(true);
            }
            pv = network.get_probability_variable(variable);
            if (pv == null) {
                System.out.println("\tInvalid variable name!");
                return(true);
            }
            if (value == null)
                pv.set_invalid_observed_index();
            else
                pv.set_observed_value(value);
            break;
        case 't':
            if (variable == null) {
                System.out.println("\tNull variable name!");
                return(true);
            }
            pv = network.get_probability_variable(variable);
            if (pv == null) {
                System.out.println("\tInvalid variable name!");
                return(true);
            }
            pv.set_explanation_index(0);
            break;
        case 'u':
            if (variable == null) {
                System.out.println("\tNull variable name!");
                return(true);
            }
            pv = network.get_probability_variable(variable);
            if (pv == null) {
                System.out.println("\tInvalid variable name!");
                return(true);
            }
            pv.set_explanation_index(-1);
            break;
        case 'i':
            if (variable == null) {
                System.out.println("\tNull variable name!");
                return(true);
            }
            pv = network.get_probability_variable(variable);
            if (pv == null) {
                System.out.println("\tInvalid variable name!");
                return(true);
            }
            inf = new Inference(network, false);
            inf.inference(variable);
            result = inf.get_result();
            System.out.print("\tPosterior marginal for " + variable + ":\n\t\t");
            for (i=0; i<result.number_values(); i++)
                System.out.print(result.get_value(i) + " ");
            System.out.println();
            break;
        case 'm':
            if (variable == null) {
                System.out.println("\tNull variable name!");
                return(true);
            }
            pv = network.get_probability_variable(variable);
            if (pv == null) {
                System.out.println("\tInvalid variable name!");
                return(true);
            }
            exp = new Expectation(network, false);
            exp.expectation(variable);
            expectation = exp.get_result();
            System.out.println("\tPosterior expectation for " + variable + ": " + expectation);
            break;
        case 'e':
            expl = new Explanation(network);
            expl.explanation();
            explanation = expl.get_explanation();
            if (explanation == null)
                System.out.println("\tNo explanation was produced!");
            else {
                for (int ie=0; ie<explanation.length; ie++)
                    if (explanation[ie][0] != null)
                        System.out.println("\tVariable " +
                            explanation[ie][0] + ": " + explanation[ie][1]);
            }
            break;
        case 'x':
            expl = new Explanation(network);
            expl.full_explanation();
            explanation = expl.get_explanation();
            if (explanation == null)
                System.out.println("\tNo explanation was produced!");
            else {
                for (int ie=0; ie<explanation.length; ie++)
                    if (explanation[ie][0] != null)
                        System.out.println("\tVariable " +
                            explanation[ie][0] + ": " + explanation[ie][1]);
            }
    }
    return(true);
  }

  /**
   * Basic help message.
   */
  static void help_message() {
    System.out.println("EBayes");
    System.out.println("Commands: ");
    System.out.println("\th -> This help message.");
    System.out.println("\tl \"name\" -> Load a network.");
    System.out.println("\to \"variable\" -> Set variable as not observed.");
    System.out.println("\to \"variable\" \"value\" -> Observe variable.");
    System.out.println("\tt \"variable\" -> Set variable as explanatory.");
    System.out.println("\tu \"variable\" -> Set variable as non-explanatory.");
    System.out.println("\ti \"variable\" -> Posterior marginal for variable.");
    System.out.println("\tm \"variable\" -> Expected value for variable.");
    System.out.println("\te -> Maximum a posteriori for explanatory variables.");
    System.out.println("\tx -> Maximum a posteriori for all variables.");
    System.out.println("\tq -> Quit.");
  }

}
