#![allow(unused)]
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

// Enum to hold all tool types
pub enum ToolVariant {
    Sassy(SassyTool),
    Swo(Swofinder),
    Chop(Chopoff),
}

pub enum Strand {
    Fwd,
    Rc,
}
pub struct Match {
    start: usize,
    end: usize,
    strand: Strand,
}

// Object safe-trait
pub trait Tool {
    fn name(&self) -> &str;

    /// Runs the tool with whatever arguments are needed.
    /// Returns the duration of just the command execution
    fn run(
        &self,
        dist: usize,
        guide_seq: &str,
        target_file_path: &str,
        out_path: &str,
        threads: usize,
    ) -> Result<Duration, String>;

    /// Parses the output by returning a vector of Match objects
    fn parse_output(&self, out_path: &str) -> Result<Vec<Match>, String>;
}

impl Tool for ToolVariant {
    fn name(&self) -> &str {
        match self {
            ToolVariant::Sassy(t) => t.name(),
            ToolVariant::Swo(t) => t.name(),
            ToolVariant::Chop(t) => t.name(),
        }
    }

    fn run(
        &self,
        dist: usize,
        guide_seq: &str,
        target_file_path: &str,
        out_path: &str,
        threads: usize,
    ) -> Result<Duration, String> {
        match self {
            ToolVariant::Sassy(t) => t.run(dist, guide_seq, target_file_path, out_path, threads),
            ToolVariant::Swo(t) => t.run(dist, guide_seq, target_file_path, out_path, threads),
            ToolVariant::Chop(t) => t.run(dist, guide_seq, target_file_path, out_path, threads),
        }
    }

    fn parse_output(&self, out_path: &str) -> Result<Vec<Match>, String> {
        match self {
            ToolVariant::Sassy(t) => t.parse_output(out_path),
            ToolVariant::Swo(t) => t.parse_output(out_path),
            ToolVariant::Chop(t) => t.parse_output(out_path),
        }
    }
}

/*
    SassyTool - just to not be confused with "sassy" code
*/

pub struct SassyTool {
    exec_path: String,
}

impl SassyTool {
    pub fn new(exec_path: &str) -> Self {
        Self {
            exec_path: exec_path.to_string(),
        }
    }
}

impl Tool for SassyTool {
    fn name(&self) -> &str {
        "SassyTool"
    }

    fn run(
        &self,
        dist: usize,
        guide_file_path: &str,
        target_file_path: &str,
        out_path: &str,
        threads: usize,
    ) -> Result<Duration, String> {
        let args = vec![
            "crispr".to_string(),
            "-g".to_string(),
            guide_file_path.to_string(),
            "-k".to_string(),
            dist.to_string(),
            "-t".to_string(),
            target_file_path.to_string(),
            "--rc".to_string(),
            "-o".to_string(),
            out_path.to_string(),
            //"--exact-suffix".to_string(),
            //"3".to_string(),
            "-j".to_string(),
            threads.to_string(),
            "-n".to_string(),
            "0.0".to_string(),
        ];

        // Execute the command and time just the execution
        let start = Instant::now();
        let mut child = Command::new(&self.exec_path)
            .args(&args)
            .stdout(Stdio::piped()) // Stream output
            .stderr(Stdio::piped()) // Stream errors
            .spawn()
            .map_err(|e| format!("Failed to spawn SassyTool process: {}", e))?;

        // Process stdout and stderr in real-time if needed
        if let Some(stdout) = child.stdout.take() {
            let stdout_reader = BufReader::new(stdout);
            for line in stdout_reader.lines() {
                if let Ok(line) = line {
                    println!("SassyTool stdout: {}", line);
                }
            }
        }

        if let Some(stderr) = child.stderr.take() {
            let stderr_reader = BufReader::new(stderr);
            for line in stderr_reader.lines() {
                if let Ok(line) = line {
                    eprintln!("SassyTool stderr: {}", line);
                }
            }
        }

        let status = child
            .wait()
            .map_err(|e| format!("Failed to wait on SassyTool process: {}", e))?;
        let duration = start.elapsed();

        if status.success() {
            Ok(duration)
        } else {
            Err(format!(
                "SassyTool command failed with exit code: {:?}",
                status.code()
            ))
        }
    }

    fn parse_output(&self, _out_path: &str) -> Result<Vec<Match>, String> {
        Ok(vec![])
    }
}

/*
    Swofinder
*/

pub struct Swofinder {
    exec_path: String,
}

impl Swofinder {
    pub fn new(exec_path: &str) -> Self {
        Self {
            exec_path: exec_path.to_string(),
        }
    }
}

impl Tool for Swofinder {
    fn name(&self) -> &str {
        "Swofinder"
    }

    fn run(
        &self,
        dist: usize,
        guide_file_path: &str,
        target_file_path: &str,
        out_path: &str,
        threads: usize,
    ) -> Result<Duration, String> {
        // Use the directory containing the Java classes as the working directory
        let working_dir = std::path::PathBuf::from(self.exec_path.clone());

        // Copy the guide file to sgRNAs.txt in the working directory
        let dest_sg_rnas_path = working_dir.join("sgRNAs.txt");
        std::fs::copy(guide_file_path, &dest_sg_rnas_path)
            .map_err(|e| format!("Failed to copy guide file: {}", e))?;

        let args = vec![
            "-cp".to_string(),
            "bin".to_string(),
            "SmithWatermanOffTarget.SmithWatermanOffTargetSearchAlign".to_string(),
            target_file_path.to_string(),
            "sgRNAs.txt".to_string(),
            out_path.to_string(),
            dist.to_string(),
            dist.to_string(),
            dist.to_string(),
            dist.to_string(),
            dist.to_string(),
            "false".to_string(),
            "0".to_string(),
            "NGG".to_string(),
            "false".to_string(),
        ];

        // Execute the command and time just the execution
        let start = Instant::now();
        let mut child = Command::new("java")
            .args(&args)
            .current_dir(working_dir)
            .stdout(Stdio::piped()) // Changed from null to piped
            .stderr(Stdio::piped()) // Changed from null to piped
            .spawn()
            .map_err(|e| format!("Failed to spawn Swofinder process: {}", e))?;

        // Stream stdout
        if let Some(stdout) = child.stdout.take() {
            let reader = BufReader::new(stdout);
            for line in reader.lines() {
                if let Ok(line) = line {
                    println!("Swofinder stdout: {}", line);
                }
            }
        }

        // Stream stderr
        if let Some(stderr) = child.stderr.take() {
            let reader = BufReader::new(stderr);
            for line in reader.lines() {
                if let Ok(line) = line {
                    eprintln!("Swofinder stderr: {}", line);
                }
            }
        }

        let status = child
            .wait()
            .map_err(|e| format!("Failed to wait on Swofinder process: {}", e))?;
        let duration = start.elapsed();

        // Clean up temporary file
        if let Err(e) = std::fs::remove_file(dest_sg_rnas_path) {
            eprintln!("Warning: Failed to remove temporary file sgRNAs.txt: {}", e);
        }

        if status.success() {
            Ok(duration)
        } else {
            Err(format!(
                "Swofinder command failed with exit code: {:?}",
                status.code()
            ))
        }
    }

    fn parse_output(&self, out_path: &str) -> Result<Vec<Match>, String> {
        Ok(vec![])
    }
}

/*
    Chopoff
*/

pub struct Chopoff {
    exec_path: String,
}

impl Chopoff {
    pub fn new(exec_path: &str) -> Self {
        Self {
            exec_path: exec_path.to_string(),
        }
    }

    // New build method specific to Chopoff
    pub fn build(
        &self,
        target_file_path: &str,
        out_dir: &str,
        distance: usize,
        db_name: &str,
        threads: usize,
    ) -> Result<Duration, String> {
        // Set JULIA_NUM_THREADS environment variable
        unsafe { std::env::set_var("JULIA_NUM_THREADS", threads.to_string()) };

        let args = vec![
            "build".to_string(),
            "--name".to_string(),
            format!("Cas9_{}", db_name),
            "--genome".to_string(),
            target_file_path.to_string(),
            "-o".to_string(),
            out_dir.to_string(),
            "--distance".to_string(),
            distance.to_string(),
            "--motif".to_string(),
            "Cas9".to_string(),
            "prefixHashDB".to_string(),
        ];

        // Execute the command and time just the execution
        let start = Instant::now();
        let mut child = Command::new(&self.exec_path)
            .args(&args)
            .stdout(Stdio::piped()) // Changed from null to piped
            .stderr(Stdio::piped()) // Changed from null to piped
            .spawn()
            .map_err(|e| format!("Failed to spawn Chopoff build process: {}", e))?;

        // Stream stdout
        if let Some(stdout) = child.stdout.take() {
            let reader = BufReader::new(stdout);
            for line in reader.lines() {
                if let Ok(line) = line {
                    println!("Chopoff build stdout: {}", line);
                }
            }
        }

        // Stream stderr
        if let Some(stderr) = child.stderr.take() {
            let reader = BufReader::new(stderr);
            for line in reader.lines() {
                if let Ok(line) = line {
                    eprintln!("Chopoff build stderr: {}", line);
                }
            }
        }

        let status = child
            .wait()
            .map_err(|e| format!("Failed to wait on Chopoff build process: {}", e))?;
        let duration = start.elapsed();

        if status.success() {
            Ok(duration)
        } else {
            Err(format!(
                "Chopoff build command failed with exit code: {:?}",
                status.code()
            ))
        }
    }
}

impl Tool for Chopoff {
    fn name(&self) -> &str {
        "Chopoff"
    }

    fn run(
        &self,
        dist: usize,
        guide_file_path: &str,
        target_file_path: &str,
        out_path: &str,
        threads: usize,
    ) -> Result<Duration, String> {
        // Create temporary guides file
        let tmp_guides_path = "tmp_guides.txt";
        let mut output_file = File::create(tmp_guides_path)
            .map_err(|e| format!("Failed to create guides file: {}", e))?;

        // Read the guide file line by line
        let guide_file =
            File::open(guide_file_path).map_err(|e| format!("Failed to open guide file: {}", e))?;
        let reader = BufReader::new(guide_file);

        // Process each line
        for line in reader.lines() {
            let line = line.map_err(|e| format!("Failed to read line: {}", e))?;
            if line.len() >= 3 {
                // Strip last 3 characters and write to output file
                let stripped = &line[..line.len() - 3];
                writeln!(output_file, "{}", stripped)
                    .map_err(|e| format!("Failed to write to guides file: {}", e))?;
            }
        }

        // Set JULIA_NUM_THREADS environment variable
        unsafe { std::env::set_var("JULIA_NUM_THREADS", threads.to_string()) };

        let args = vec![
            "search".to_string(),
            "--database".to_string(),
            target_file_path.to_string(),
            "--guides".to_string(),
            tmp_guides_path.to_string(),
            "--output".to_string(),
            out_path.to_string(),
            "--distance".to_string(),
            dist.to_string(),
            "prefixHashDB".to_string(),
        ];

        // Execute the command and stream output
        let start = Instant::now();
        let mut child = Command::new(&self.exec_path)
            .args(&args)
            .stdout(Stdio::piped()) // Pipe stdout
            .stderr(Stdio::piped()) // Pipe stderr
            .spawn()
            .map_err(|e| format!("Failed to spawn Chopoff search process: {}", e))?;

        // Stream stdout
        if let Some(stdout) = child.stdout.take() {
            let reader = BufReader::new(stdout);
            for line in reader.lines() {
                if let Ok(line) = line {
                    println!("Chopoff stdout: {}", line);
                }
            }
        }

        // Stream stderr
        if let Some(stderr) = child.stderr.take() {
            let reader = BufReader::new(stderr);
            for line in reader.lines() {
                if let Ok(line) = line {
                    eprintln!("Chopoff stderr: {}", line);
                }
            }
        }

        let status = child
            .wait()
            .map_err(|e| format!("Failed to wait on Chopoff search process: {}", e))?;
        let duration = start.elapsed();

        // Clean up temporary file
        if let Err(e) = std::fs::remove_file(tmp_guides_path) {
            eprintln!(
                "Warning: Failed to remove temporary file {}: {}",
                tmp_guides_path, e
            );
        }

        if status.success() {
            Ok(duration)
        } else {
            Err(format!(
                "Chopoff search command failed with exit code: {:?}",
                status.code()
            ))
        }
    }

    fn parse_output(&self, out_path: &str) -> Result<Vec<Match>, String> {
        // Simplified for this example
        Ok(vec![])
    }
}
