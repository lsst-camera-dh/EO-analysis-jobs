#!/usr/bin/env python
"""
Validator script for BOT bright and dark pixel analyses.
"""
from bot_eo_validators import run_validator
run_validator('bright_defects')
run_validator('dark_defects')