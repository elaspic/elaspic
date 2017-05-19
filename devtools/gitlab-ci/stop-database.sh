#!/bin/bash

set -ev

ps aux | grep mysql | grep -v grep | awk '{print $2}' | xargs -i{} kill {}
